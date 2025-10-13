
## LDpred function per chromosome
LDpredest <-
  function(df_gwas, df_map, LD_path="~/Downloads/LDpred/ldref_hm3_plus/",
          prefix = "LD_with_blocks_chr", MAF_threshold = NA,
          h2input="", NCORES=1){
    # colnames(df_gwas) = c("rsid","chr","pos","beta","se","N")
    # colnames(df_map) = c("rsid","chr","pos","ld","MAF")
    ## dependencies
    require(bigsnpr)
    require(dplyr)
    require(Matrix)

    time_start <- date()
    cat(paste0("Analysis started on ",time_start,"\n"))
    
    ## estimate LDSC heritability if missing "h2input"
    if (h2input == "") {
      st = Sys.time()
      cat(paste0("Estimating LDSC heritability ...\t"))
      df = df_gwas %>% inner_join(df_map)
      LDscore = df$ld
      res_ldsc <- with(df, bigsnpr::snp_ldsc(LDscore, length(LDscore), chi2 = (beta/se)^2, sample_size = N, blocks = NULL))
      h2input = res_ldsc[["h2"]]
      ed = Sys.time()
      cat(paste0("Time ",difftime(ed,st,units="mins")," min \n"))
      cat(paste0("LDSC estimates of heritability = ",round(h2input,8),"\n"))
    }
    if(is.na(h2input) | h2input<0 | h2input>1){
      error.message <- "The input of h2input has to be a number between 0 and 1. \n"
      stop(error.message)
    }
    ## SNP info in LD
    if (is.numeric(MAF_threshold) && length(MAF_threshold) == 1) {
      df_map <- df_map %>% select(rsid,chr,pos,MAF)
      } else {
        df_map <- df_map %>% select(rsid,chr,pos)
      }
    ## GWAS
    df_gwas <- df_gwas %>%
      select(rsid,chr,pos,beta,se,N) %>%
      mutate(beta = ifelse(is.na(beta),0,beta), 
             se = ifelse(is.na(se),1,se), 
             N = ifelse(is.na(N),median(df_gwas$N,na.rm = TRUE),N), 
             scale = sqrt(beta^2+N*se^2), 
             sbeta = beta/scale)

    ## function: estimation per chromosome
    estimation_chr <- function(df_gwas_chr, df_map_chr, R_chr, h2input, M, NCORES){
      # cat(paste0("estimation winID ", j, "/", ceiling(nrow(df_gwas_chr)/ws)))
      R_chr <- R_chr[df_map_chr$mapID, df_map_chr$mapID]
      R_SFBM <- bigsparser::as_SFBM(R_chr)
      res_chr <- list(h2input = h2input, M = M, M_chr = ncol(R_chr))
      
      ## LDpred2-inf
      res_chr$b_inf_LDpred <- bigsnpr::snp_ldpred2_inf(corr = R_SFBM, df_beta = df_gwas_chr %>% rename("beta_se"="se","n_eff"="N"), 
                                                       h2 = h2input*nrow(df_gwas_chr)/M) / df_gwas_chr$scale
      res_chr$h2_inf_LDpred <- as.numeric(res_chr$b_inf_LDpred %*% R_chr %*% res_chr$b_inf_LDpred)
      res_chr$Rb_inf_LDpred <- as.numeric(R_chr %*% res_chr$b_inf_LDpred)
      res_chr$h2_inf_LDpred <- as.numeric(res_chr$b_inf_LDpred %*% R_chr %*% res_chr$b_inf_LDpred)
      ## LDpred2-auto
      model <- bigsnpr::snp_ldpred2_auto(corr = R_SFBM, df_beta = df_gwas_chr %>% rename("beta_se"="se","n_eff"="N"), 
                                         h2_init = h2input*nrow(df_gwas_chr)/M, ncores = NCORES, 
                                         vec_p_init = bigsnpr::seq_log(1e-4,0.2,length.out=30), 
                                         use_MLE= FALSE,  ### this is for convergence issues or GWAS power is low 
                                         allow_jump_sign = FALSE, shrink_corr = 0.95)
      range <- sapply(model, function(auto) diff(range(auto$corr_est)))
      keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
      
      cat("chr", unique(df_gwas_chr$chr), " kept_chains:", length(keep), "\n") ###### Sanity CHECK  
      
      # 👇 Diagnostics
      cat("chr", unique(df_gwas_chr$chr), 
          " kept_chains:", length(keep), 
          " range summary:", paste0(round(summary(range), 3), collapse = " "), "\n")
      
      # Quick look at model object
      cat("chr", unique(df_gwas_chr$chr), 
          " model length:", length(model), 
          " beta_est length (first chain):", length(model[[1]]$beta_est), "\n")
      
      res_chr$b_auto_LDpred <- rowMeans(sapply(model[keep], function(auto) auto$beta_est)) / df_gwas_chr$scale
      res_chr$pp_auto_LDpred <- rowMeans(sapply(model[keep], function(auto) auto$postp_est))
      res_chr$pi_auto_LDpred <- mean(sapply(model[keep], function(auto) auto$p_est))
      res_chr$h2_auto_LDpred <- as.numeric(res_chr$b_auto_LDpred %*% R_chr %*% res_chr$b_auto_LDpred)
      res_chr$Rb_auto_LDpred <- as.numeric(R_chr %*% res_chr$b_auto_LDpred)
      res_chr$h2_auto_LDpred <- as.numeric(res_chr$b_auto_LDpred %*% R_chr %*% res_chr$b_auto_LDpred)
      rm(model,range,keep)
      
      return(as.data.frame(res_chr))
    }
    ## estimate genetic effects
    st = Sys.time()
    cat(paste0("Estimating genetic effects ...\n"))
    df_est <- NULL
    for (CHR in sort(unique(intersect(df_gwas$chr, df_map$chr)))) {
      st_chr = Sys.time()
      cat(paste0("Estimating in chromosome ",CHR," ...\t"))
      ## data
      df_gwas_chr <- df_gwas %>% filter(chr==CHR)
      df_map_chr <- df_map %>% filter(chr==CHR) %>% 
        mutate(mapID = rank(pos)) %>% 
        distinct(rsid, .keep_all = TRUE)
      if (is.na(MAF_threshold)) {
        df_map_chr <- inner_join(df_map_chr,df_gwas_chr, by = c("rsid","chr","pos"))
      } else if (is.numeric(MAF_threshold) && length(MAF_threshold) == 1) {
        cat("removing SNPs with MAF <",MAF_threshold,"... \t")
        df_map_chr <- filter(df_map_chr,MAF >= MAF_threshold) %>% 
          inner_join(df_gwas_chr, by = c("rsid","chr","pos")) %>% 
          select(-MAF) 
        df_gwas_chr <- df_gwas_chr %>% 
          inner_join(df_map_chr %>% select(rsid,chr,pos), by = c("rsid","chr","pos"))
        } else{
          stop("MAF_threshold should be numeric.")
        }
      cat("chr", CHR, "nSNPs after filtering:", nrow(df_gwas_chr), "\n") ####### SANITY CHECK 
      R_chr <- readRDS(paste0(LD_path,prefix,CHR,".rds"))
      R_chr@x[is.na(R_chr@x)] <- 1e-10
      ## computation
      df_est_chr <- estimation_chr(df_gwas_chr = df_gwas_chr, df_map_chr = df_map_chr, 
                                   R_chr = R_chr, h2input = h2input, M = nrow(df_gwas), 
                                   NCORES = NCORES)
      df_est_chr <- cbind(df_gwas_chr, df_est_chr)
      df_est <- rbind(df_est, df_est_chr)
      rm(df_gwas_chr, df_map_chr, R_chr, df_est_chr)
      ed_chr = Sys.time()
      cat(paste0("Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
    }
    ed = Sys.time()
    
    cat(paste0("LDpred estimation is DONE with Time ",difftime(ed,st,units="mins")," min \n"))
    
    cat(paste0("\n"))
    
    time_end <- date()
    cat(paste0("Analysis finished on ",time_end,"\n"))
    
    return(df_est)
  }


