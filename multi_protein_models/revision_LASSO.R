library(glmnet)
library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(survival)
library(ggplot2)
library(patchwork)
library(pROC)
library(caret)
library(arrow)
library(doParallel)
library(foreach)
library(stats)
library(stabs, lib.loc="~/R_library/")
library(predtools, lib = "~/R_library/")
library(ggrepel)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

## Run stability selection when selecting features with LASSO - logistic regression and AUC extraction  - And calibration 
lasso_with_stability <- function(X, Y, seed = 1234, testsize = 0.3,
                                 cutoff = 0.6, PFER = 3, B = 100,
                                 nfolds_cv = 5, nlambda = 200,
                                 ece_bins = 15     # number of bins for ECE
) {
  
  # ---- output containers ----
  results <- data.frame(
    Disease   = colnames(Y),
    AUC       = NA_real_,
    AUC_Lower = NA_real_,
    AUC_Upper = NA_real_,
    n_features  = NA_integer_,
    features    = NA_character_,
    # calibration / proper scoring
    ECE        = NA_real_,
    stringsAsFactors = FALSE
  )
  
  pred_list <- vector("list", ncol(Y))
  
  # ---- parallel backend ----
  cl <- makeCluster(parallelly::availableCores(omit = 1))
  registerDoParallel(cl)
  
  out <- foreach(i = 1:ncol(Y),
                 .packages = c("caret","glmnet","pROC","stats")) %dopar% {
                   
                   
                   # =========================
                   # Calibration metrics 
                   # =========================
                   
                   # ECE (Expected Calibration Error) via binning:
                   # sum_k (n_k/N) * | mean(y)_k - mean(p)_k |
                   ece_binary <- function(y, p, n_bins = 15) {
                     y <- as.numeric(y); p <- as.numeric(p)
                     ok <- is.finite(y) & is.finite(p)
                     y <- y[ok]; p <- p[ok]
                     if (length(p) == 0) return(NA_real_)
                     
                     p <- pmin(pmax(p, 0), 1)
                     bins <- seq(0, 1, length.out = n_bins + 1)
                     ece <- 0
                     N <- length(p)
                     
                     for (i in seq_len(n_bins)) {
                       # include right edge in last bin
                       msk <- if (i < n_bins) {
                         (p > bins[i]) & (p <= bins[i + 1])
                       } else {
                         (p > bins[i]) & (p <= bins[i + 1] + 1e-9)
                       }
                       if (!any(msk)) next
                       acc_bin  <- mean(y[msk] == 1)
                       conf_bin <- mean(p[msk])
                       ece <- ece + (sum(msk) / N) * abs(acc_bin - conf_bin)
                     }
                     as.numeric(ece)
                   }
                   
                   
                   .libPaths("~/R_library/")
                   library(stabs)
                   
                   set.seed(seed)
                   
                   
                   # ---- outcome ----
                   y0 <- Y[, i] ## extract the i-th outcome
                   y0 <- y0[!is.na(y0)]  # drop NA ( which are dropped indv but present in other diseases) 
                   
                   if (length(unique(y0)) <= 1L) return(NULL)
                   
                   # align X to y
                   X_filt <- X[names(y0), , drop = FALSE]
                   
                   # safety: ensure binary 0/1 numeric
                   y <- as.integer(as.character(y0))
                   if (any(!y %in% c(0L, 1L))) {
                     # if y0 is logical, factor, etc.
                     y <- as.integer(y0)
                   }
                   if (any(!y %in% c(0L, 1L))) return(NULL)
                   
                   # ---- train/test split ----
                   train_idx <- caret::createDataPartition(y, p = 1 - testsize, list = FALSE)
                   
                   y_train <- y[train_idx] ### endpoint splitting
                   y_test  <- y[-train_idx]
                   
                   X_train <- X_filt[train_idx, , drop = FALSE] ### proteins splitting
                   X_test  <- X_filt[-train_idx, , drop = FALSE]
                   
                   
                   # ---- stability selection on TRAIN only ----
                   set.seed(seed + i)
                   
                   stab_fit <- stabs::stabsel(
                     x = X_train,
                     y = y_train,
                     fitfun = stabs::glmnet.lasso,
                     args.fitfun = list(family = "binomial", alpha = 1),
                     cutoff = cutoff,
                     PFER   = PFER,
                     B      = B 
                   )
                   
                   stable_idx <- stab_fit$selected
                   
                   stable_vars <- if (length(stable_idx)) colnames(X_train)[stable_idx] else character(0)
                   
                   # if nothing stable, return NA metrics + empty feature list
                   if (length(stable_idx) == 0L) {
                     return(list(
                       i = i,
                       auc = c(NA_real_, NA_real_, NA_real_),
                       ECE=NA_real_,
                       final_vars = character(0),
                       y_test = y_test,
                       pred = rep(NA_real_, length(y_test)),
                       stable_vars = stable_vars,
                       selprob = stab_fit$max,
                       final_weights = numeric(0)
                     ))
                   }
                   
                   # ---- SPECIAL CASE: only 1 stable feature -> glmnet will error ----
                   if (length(stable_idx) == 1L) {
                     
                     x1_train <- as.numeric(X_train[, stable_idx])
                     x1_test  <- as.numeric(X_test[,  stable_idx])
                     
                     # estimate weight for this single feature on TRAIN (unpenalized)
                     fit1 <- glm(y_train ~ x1_train, family = "binomial")
                     w_final <- as.numeric(coef(fit1)[2])          # slope for the feature
                     
                     # build score (single-feature score)
                     score_train <- x1_train * w_final
                     score_test  <- x1_test  * w_final
                     
                     # final risk model 
                     model <- glm(y_train ~ score_train, family = "binomial")
                     pred <- predict(model,
                                     newdata = data.frame(score_train = score_test),
                                     type = "response")
                     
                     # discrimination
                     auc_val <- as.numeric(pROC::auc(y_test, pred))
                     ci_val  <- as.numeric(pROC::ci.auc(y_test, pred))
                     
                     # calibration metrics 
                     ece   <- ece_binary(y_test, pred, n_bins = ece_bins)
                     
                     return(list(
                       i = i,
                       auc = c(auc_val, ci_val[1], ci_val[3]),
                       ECE = ece,
                       final_vars = colnames(X_train)[stable_idx],
                       y_test = y_test,
                       pred = pred,
                       stable_vars = stable_vars,
                       selprob = stab_fit$max,
                       final_weights = w_final
                     ))
                   }
                   
                   # ---- refit final model on TRAIN using only stable vars ----
                   cv_fit2 <- glmnet::cv.glmnet(
                     X_train[, stable_idx, drop = FALSE],
                     y_train,
                     family = "binomial",
                     alpha = 1,
                     nfolds = nfolds_cv,
                     nlambda = nlambda,
                     type.measure = "auc"
                   )
                   
                   # take lambda.min
                   coefs2 <- coef(cv_fit2, s = "lambda.min")
                   w2 <- as.numeric(coefs2[-1, , drop = FALSE])  # weights for stable vars
                   sel2 <- which(w2 != 0)
                   
                   if (length(sel2) == 0L) {
                     # If lambda.min shrinks all to 0, try lambda.1se as a conservative alternative
                     coefs2 <- coef(cv_fit2, s = "lambda.1se")
                     w2 <- as.numeric(coefs2[-1, , drop = FALSE])
                     sel2 <- which(w2 != 0)
                   }
                   
                   # if still none selected, return NA
                   if (length(sel2) == 0L) {
                     return(list(
                       i = i,
                       auc = c(NA_real_, NA_real_, NA_real_),
                       final_vars = character(0),
                       pred = rep(NA_real_, length(y_test)),
                       y_test = y_test,
                       stable_vars = stable_vars,
                       selprob = stab_fit$max
                     ))
                   }
                   
                   
                   # map sel2 (within stable set) back to original column indices
                   cols_final <- stable_idx[sel2]
                   w_final <- w2[sel2]
                   final_vars <- colnames(X_train)[cols_final]
                   
                   #---- build linear score 
                   score_train <- as.numeric(X_train[, cols_final, drop = FALSE] %*% w_final)
                   score_test  <- as.numeric(X_test[, cols_final, drop = FALSE] %*% w_final)
                   
                   score = score_train
                   
                   # ---- logistic regression with score on training set----
                   model <- glm(y_train ~ score, family = "binomial")
                   
                   
                   # prediction on test set
                   pred <- predict(model,
                                   newdata = data.frame(score = score_test),
                                   type = "response")
                   
                   # ---- AUC + CI ----
                   auc_val <- as.numeric(pROC::auc(y_test, pred))
                   ci_val  <- as.numeric(pROC::ci.auc(y_test, pred))
                   
                   # -------- calibration on TEST --------
                   ece   <- ece_binary(y_test, pred, n_bins = ece_bins)
                   
                   list(
                     i = i,
                     auc = c(auc_val, ci_val[1], ci_val[3]),
                     ECE=ece,
                     final_vars = final_vars,
                     pred = pred,
                     y_test = y_test,
                     stable_vars = stable_vars,
                     selprob = stab_fit$max,
                     final_weights = w_final
                   )
                 }
  
  stopCluster(cl)
  
  # ---- collect ----
  for (piece in out) {
    if (is.null(piece)) next
    i <- piece$i
    
    results$AUC[i]       <- piece$auc[1]
    results$AUC_Lower[i] <- piece$auc[2]
    results$AUC_Upper[i] <- piece$auc[3]
    
    results$ECE[i]   <- piece$ECE
    
    if (length(piece$final_vars) > 0) {
      results$n_features[i] <- length(piece$final_vars)
      results$features[i]   <- paste(piece$final_vars, collapse = ",")
    } else {
      results$n_features[i] <- 0L
      results$features[i]   <- ""
    }
    
    pred_list[[i]] <- piece
  }
  
  list(results = results, predictions = pred_list)
}


#### LOAD DATA ####

### file with first occurrences of INTERVENE phenotypes 
intervene_fi= arrow::read_parquet("/scratch/project_2007428/projects/prj_010_phenotype_file/share/ukb78537_INTERVENE_phenotype-2024-08-01.parquet")
intervene_fi= intervene_fi %>% select(ID,SEX,END_OF_FOLLOWUP, contains("_DATE"))

### filter phenotypes for the intervene. note that from the 39 endpoints, BMI and Covid-19 severity are not present as first occurrences, so the final number is 37 
definitions= fread("data/Intervene_flagship_endpoint_collection_Definitions.tsv") 
definitions = definitions %>% select(Endpoint,`FinnGen endpoint`)

## load age and sex covariates, first assessment visit date and death
phenos=fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv")
selphenos= phenos %>% select(ID=eid, age='21022-0.0',
                             sex='31-0.0', baseline='53-0.0', death='40000-0.0')
selphenos = selphenos %>% filter(ID %in% intervene_fi$ID)

intervene_fi= merge(intervene_fi, selphenos, by="ID")

allpheno= fread("data/idefix_olinksoma/pheno_residuals.tsv")

intervene_fi= intervene_fi %>% filter(ID %in% allpheno$eid)

### creating a matrix with all the endpoints 
allendpoints = data.frame(ID = intervene_fi$ID)
for (dis in colnames(intervene_fi %>% select(contains("_DATE")))) {
  df= intervene_fi %>% select(ID, END_OF_FOLLOWUP, dis, baseline) 
  df_merge= df %>% mutate(fol = difftime(!!as.name(dis), baseline, units = "weeks")/52.25) ###note that I converted in years!!
  
  ## excluding all prevalent cases or incident cases recorded within the first 6 months of follow-up
  df_filt= df_merge %>% filter(fol>0.5 | is.na(fol))
  ### binary outcome, 1= incidence, 0 = right censored
  df_bin= df_filt %>% mutate(!!as.name(dis) := ifelse(is.na(!!as.name(dis)),0,1)) %>% select(ID, all_of(dis))
  
  allendpoints <- full_join(allendpoints, df_bin, by = "ID")
}
rownames(allendpoints)= allendpoints$ID
allendpoints$ID= NULL
m_allendpoints= as.matrix(allendpoints)

### scale, and create matrix
load_scale_filter_matrix <- function(file, match_ids, id_col = 1) {
  dt <- fread(file)
  setDT(dt)[, (names(dt)[-id_col]) := lapply(.SD, scale), .SDcols = names(dt)[-id_col]]
  
  mat <- as.matrix(dt, rownames = TRUE)
  mat <- mat[rownames(mat) %in% match_ids, ]
  return(mat)
}

#------------------------------------------------ LOAD ALL THE PREDICTORS ---------------------------------------
m_allpheno= load_scale_filter_matrix("data/idefix_olinksoma/pheno_residuals.tsv", rownames(m_allendpoints))
m_resid= load_scale_filter_matrix("data/idefix_olinksoma/residuals_output_nocov.tsv", rownames(m_allendpoints))
# Drop column after loading
m_resid= m_resid[, colnames(m_resid) != "AggAbsResid"]

m_cisresid= load_scale_filter_matrix("data/cis_residuals_output_nocov.tsv", rownames(m_allendpoints))
m_cisresid= m_cisresid[, colnames(m_cisresid) != "AggAbsResid"]

m_filtcisresiduals= load_scale_filter_matrix("data/filtered_cis_pPGS/cis_residuals_output_nocov.tsv", rownames(m_allendpoints))
m_filtcisresiduals= m_filtcisresiduals[, colnames(m_filtcisresiduals) != "AggAbsResid"]

m_filtresiduals= load_scale_filter_matrix("data/filtered_gw_pPGS/gw_residuals_output_nocov.tsv", rownames(m_allendpoints))
m_filtresiduals= m_filtresiduals[, colnames(m_filtresiduals) != "AggAbsResid"]

##------------------------------ RUN LASSO ----------------------------------------------------------------------------

#################### run LASSO   

### proteins
pred_prot= lasso_with_stability(m_allpheno,m_allendpoints)

### residuals
pred_resid= lasso_with_stability(m_resid, m_allendpoints)

##### cis-residuals
pred_cisresid= lasso_with_stability(m_cisresid,m_allendpoints)


#### create a matrix with filt residuals + non filt residuals ( in case there are not filtered )
## and run LASSO with the combination of filt + non filt
make_residual_matrix_for_endpoint <- function(endpoint,
                                              m_allresiduals,
                                              m_filt_residuals,
                                              base_prefix = "StandResid_",
                                              sep = "__") {
  
  X <- as.matrix(m_allresiduals)
  base_cols = colnames(m_allresiduals)  # StandResid_<protein>
  proteins  = sub(paste0("^", base_prefix), "", base_cols)
  
  filt_cols= paste0(base_prefix, proteins, sep, endpoint)  # StandResid_<protein>__I9_CHD
  
  has = filt_cols %in% colnames(m_filt_residuals)
  if (any(has)) {
    # align rows by rownames 
    common_ids <- intersect(rownames(X), rownames(m_filt_residuals))
    X <- X[common_ids, , drop = FALSE]
    X[, base_cols[has]] = as.matrix(m_filt_residuals[common_ids, filt_cols[has], drop = FALSE])
  }
  
  colnames(X) <- base_cols
  X
}

run_lasso_all_endpoints_filtered <- function(m_allendpoints,
                                             m_allresiduals,
                                             m_filt_residuals,
                                             ...) {
  endpoints <- colnames(m_allendpoints)
  
  fits <- vector("list", length(endpoints))
  names(fits) <- endpoints
  
  for (ep in endpoints) {
    X_ep <- make_residual_matrix_for_endpoint(
      endpoint = ep,
      m_allresiduals = m_allresiduals,
      m_filt_residuals = m_filt_residuals
    )
    
    Y_ep <- m_allendpoints[rownames(X_ep), ep, drop = FALSE]
    colnames(Y_ep) <- ep  
    
    # your lasso function stays exactly as it is:
    fits[[ep]] <- lasso_with_stability(X = X_ep, Y = Y_ep, ...)
  }
  
  fits
}


lasso_residuals_filtered <- run_lasso_all_endpoints_filtered(
  m_allendpoints   = m_allendpoints,
  m_allresiduals   = m_resid,
  m_filt_residuals = m_filtresiduals)


lasso_cisresiduals_filtered <- run_lasso_all_endpoints_filtered(
  m_allendpoints   = m_allendpoints,
  m_allresiduals   = m_resid,
  m_filt_residuals = m_filtcisresiduals)

#-------------------------------------
######### extract dataframes
#-------------------------------------
extract_auc_results <- function(model_result, model_name) {
  df <- as.data.frame(model_result[[1]])
  colnames(df) <- gsub("^results\\.", "", colnames(df))
  df$Model <- model_name
  
  preds <- model_result[[2]]
  pred_long <- vector("list", nrow(df)) 
  
  for (i in seq_len(nrow(df))) {
    pi <- preds[[i]]
    if (is.null(pi) || is.null(pi$y_test)) next
    
    y <- pi$y_test
    pred <- pi$pred
    
    df$n_total[i]    <- length(y)
    df$n_cases[i]    <- sum(y == 1)
    df$n_controls[i] <- sum(y == 0)
    
    ## individual-level table
    pred_long[[i]] <- data.frame(
      Disease = df$Disease[i],
      Model   = model_name,
      y       = y,
      pred    = pred
    )
  }
  
  list(
    summary = df,
    predictions = bind_rows(pred_long)
  )
}

prot_auc = extract_auc_results(pred_prot, "unadjusted proteins")

resid_auc = extract_auc_results(pred_resid, "genetically adjusted proteins")

cisresid_auc= extract_auc_results(pred_cisresid, "cis-genetically adjusted proteins")

#### extract combined dataframes '
extract_auc_results_from_endpoint_list <- function(fits_by_endpoint, model_name) {
  summary_list <- vector("list", length(fits_by_endpoint))
  preds_long   <- vector("list", length(fits_by_endpoint))
  
  eps <- names(fits_by_endpoint)
  if (is.null(eps)) eps <- seq_along(fits_by_endpoint)
  
  for (k in seq_along(fits_by_endpoint)) {
    ep_res <- fits_by_endpoint[[k]]
    if (is.null(ep_res)) next
    
    # ---- summary (1 row) ----
    df <- as.data.frame(ep_res$results, stringsAsFactors = FALSE)
    if (nrow(df) == 0) next
    df$Model <- model_name
    summary_list[[k]] <- df
    
    # ---- predictions (list length 1; may be NULL if endpoint skipped) ----
    pi <- ep_res$predictions[[1]]
    if (is.null(pi) || is.null(pi$y_test)) next
    
    y <- pi$y_test
    pred <- pi$pred
    
    df$n_total    <- length(y)
    df$n_cases    <- sum(y == 1)
    df$n_controls <- sum(y == 0)
    
    # write back counts into the same df row
    summary_list[[k]] <- df
    
    preds_long[[k]] <- data.frame(
      Disease = df$Disease[1],
      Model   = model_name,
      y       = y,
      pred    = pred
    )
  }
  
  list(
    summary      = dplyr::bind_rows(summary_list),
    predictions  = dplyr::bind_rows(preds_long)
  )
}

residFILT_auc= extract_auc_results_from_endpoint_list(lasso_residuals_filtered, model_name = "filtered-genetically adjusted proteins")
cisresidFILT_auc= extract_auc_results_from_endpoint_list(lasso_cisresiduals_filtered, model_name = "filtered-cis-genetically adjusted proteins")


all_models = list(prot_auc, resid_auc, cisresid_auc, residFILT_auc, cisresidFILT_auc)

saveRDS(all_models, "data/revised_LASSO_AUC_all_models_extracted.rds")


### All predictions
allpred= rbind(prot_auc$summary, resid_auc$summary, cisresid_auc$summary, 
               residFILT_auc$summary, cisresidFILT_auc$summary)
allpred$Disease= gsub("_DATE","",allpred$Disease)
allpred= merge(allpred,definitions,
               by.x="Disease", by.y="FinnGen endpoint")
allpred$features= gsub("StandResid_","",allpred$features)

write.table(allpred, "data/revision_LASSO_AUC_summarydf.tsv",quote=F, sep="\t",row.names = F)


allpred= fread("data/revision_LASSO_AUC_summarydf.tsv")
### standard error
allpred= allpred %>% mutate(SE= (AUC_Upper-AUC_Lower)/3.92)


# ---- choose models once ----

### data-frame , metrics and plot between 2 Models 
compare_auc_models <- function(allpred, models_keep, compare = models_keep) {
  stopifnot(length(compare) == 2)
  
  m1 <- compare[1]  # baseline
  m2 <- compare[2]  # new
  
  # ---- filter ----
  allpred_f <- allpred %>% filter(Model %in% models_keep)
  
  # ---- wide ----
  df_wide <- allpred_f %>%
    pivot_wider(
      id_cols = c(Disease, Endpoint),
      names_from = Model,
      values_from = c(AUC, AUC_Lower, AUC_Upper, SE, ECE, n_features, features),
      names_glue = "{.value}_{Model}"
    )
  
  # dynamic cols
  AUC1 <- paste0("AUC_", m1); AUC2 <- paste0("AUC_", m2)
  SE1  <- paste0("SE_",  m1); SE2  <- paste0("SE_",  m2)
  
  # ---- z-test (one-tailed: m2 > m1) ----
  df_wide <- df_wide %>%
    mutate(
      zscore = (.data[[AUC2]] - .data[[AUC1]]) / sqrt(.data[[SE2]]^2 + .data[[SE1]]^2),
      pval_z = 1 - pnorm(zscore),
      sig_zscore = if_else(pval_z < 0.05, "yes", "no"),
      diff = .data[[AUC2]] - .data[[AUC1]],
      perc_gain = (diff / .data[[AUC1]]) * 100
    )
  
  # ---- join back for plotting ----
  allpred_out <- allpred_f %>%
    left_join(
      df_wide %>% select(Disease, Endpoint, sig_zscore),
      by = c("Disease", "Endpoint")
    ) %>%
    mutate(
      AUC_Lower = pmax(AUC_Lower, 0.5),
      AUC       = pmax(AUC, 0.5),
      alpha_ztest = if_else(sig_zscore == "yes", 1, 0.5, missing = 0.5)
    )
  
  
  endpoint_levels <- allpred_out %>%
    group_by(Endpoint) %>%
    summarise(Max_AUC = max(AUC, na.rm = TRUE), .groups = "drop") %>%
    arrange(Max_AUC) %>%              # low -> high
    pull(Endpoint)
  
  allpred_out <- allpred_out %>%
    group_by(Endpoint) %>%
    mutate(Max_AUC = max(AUC, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      Endpoint = reorder(Endpoint, Max_AUC),
      Endpoint = factor(Endpoint, levels = endpoint_levels)
    )

  
  # ---- alternating background (based on displayed Endpoint order) ----
  stripe_df <- allpred_out %>%
    distinct(Endpoint) %>%
    mutate(y = as.numeric(Endpoint)) %>%      # numeric position of the factor level
    filter(y %% 2 == 1) %>%                   # stripe odd rows
    transmute(ymin = y - 0.5, ymax = y + 0.5)
  
  # manual colors
  cols <- if (length(models_keep) == 2) {
    setNames(c("#4DA1D0", "#E69F00"), models_keep)
  } else {
    setNames(scales::hue_pal()(length(models_keep)), models_keep)
  }
  
  p <- ggplot(allpred_out, aes(
    y = Endpoint, x = AUC,
    xmin = AUC_Lower, xmax = AUC_Upper,
    color = Model, alpha = alpha_ztest
  )) +
    
    geom_rect(
      data = stripe_df,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey80", alpha = 0.3
    ) +
    
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    geom_errorbar(width = 0.2, position = position_dodge(width = 0.6)) +
    scale_alpha_identity() +
    geom_text(
      aes(x = 0.5, label = ifelse(sig_zscore == "yes", "*", NA)),
      size = 7, color = "black", vjust = 0.8
    ) +
    coord_cartesian(xlim = c(0.5, 0.9)) +
    scale_color_manual(values = cols) +
    labs(x = "AUC (95% CI)", y = "Disease") +
    annotate(
      "text", y = 1.5, x = 0.7,
      label = paste0("* One-tailed p < 0.05 (ΔAUC z-score)"),
      hjust = 0, size = 3.5
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 12))
  
  #-----------------------------------------------
  ##### how many features overlap between models ? 
  #----------------------------------------------
  
  parse_features <- function(x) {
    x <- as.character(x)
    if (is.na(x) || trimws(x) == "") return(character(0))
    feats <- unlist(strsplit(x, ","))
    feats <- trimws(feats)
    feats <- feats[feats != ""]
    unique(feats)
  }
  
  build_overlap_from_feature_df <- function(df, model_unadj, model_adj) {
    df_w <- df %>%
      filter(Model %in% c(model_unadj, model_adj)) %>%
      select(Endpoint, Model, features) %>%
      distinct() %>%
      tidyr::pivot_wider(names_from = Model, values_from = features)
    
    df_w %>%
      rowwise() %>%
      mutate(
        feats_unadj   = list(parse_features(.data[[model_unadj]])),
        feats_adj     = list(parse_features(.data[[model_adj]])),
        N_shared      = length(intersect(feats_unadj, feats_adj)),
        N_unadj_only  = length(setdiff(feats_unadj, feats_adj)),
        N_adj_only    = length(setdiff(feats_adj, feats_unadj))
      ) %>%
      ungroup() %>%
      select(Endpoint, N_unadj_only, N_shared, N_adj_only)
  }
  
  overlap_df <- build_overlap_from_feature_df(
    df = allpred,          # use full df with features
    model_unadj = m1,
    model_adj   = m2
  )
  
  plot_df <- overlap_df %>%
    mutate(Endpoint = factor(Endpoint, levels = endpoint_levels)) %>%
    filter(!is.na(Endpoint)) %>%
    pivot_longer(
      cols = c(N_unadj_only, N_shared, N_adj_only),
      names_to = "Category",
      values_to = "N_proteins"
    ) %>%
    mutate(
      Category = factor(
        Category,
        levels = c("N_unadj_only", "N_shared", "N_adj_only"),
        labels = c(paste0(m1, " only"), "shared", paste0(m2, " only"))
      )
    )
  
  p_right <- ggplot(plot_df, aes(y = Endpoint, x = N_proteins, fill = Category)) +
    geom_col(width = 0.6, alpha = 0.7) +
    scale_fill_manual(
      values = setNames(
        c("blue4", cols[[m2]], cols[[m1]]),
        c("shared", paste0(m2, " only"), paste0(m1, " only"))
      )
    ) +
    scale_y_discrete(limits = endpoint_levels, expand = expansion(add = 0)) +
    scale_x_continuous(breaks = 1:11, minor_breaks = NULL) +
    labs(y = NULL, x = "Number of selected features", fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),   # kill horizontal lines
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    )
  
  list(
    df_wide = df_wide,
    allpred = allpred_out,
    mean_perc_gain_sig = df_wide %>% filter(sig_zscore == "yes") %>% pull(perc_gain) %>% mean(na.rm = TRUE),
    p_left = p,
    p_right = p_right,
    endpoint_levels = endpoint_levels
  )
}


prot_resid_models <- compare_auc_models(allpred, c("unadjusted proteins","genetically adjusted proteins"))
prot_cisresid_models <- compare_auc_models(allpred, c("unadjusted proteins","cis-genetically adjusted proteins"))
prot_filtcisresid_models <- compare_auc_models(allpred, c("unadjusted proteins","filtered-cis-genetically adjusted proteins"))

#-----------------------------------------------
##### plot LASSO AUC
#----------------------------------------------

### FIGURE 5
pdf("plots/revised_Figure5.pdf", width = 13.5, height = 9)
prot_resid_models$p_left + prot_resid_models$p_right
dev.off()

png("plots/revised_Figure5.png", width = 13.5, height = 9, units = "in", res = 300)
prot_resid_models$p_left + prot_resid_models$p_right
dev.off()


#### cis - figure 5 
pdf("plots/revised_cisLASSO.pdf", width = 13.5, height = 9)
prot_cisresid_models$p_left + prot_cisresid_models$p_right
dev.off()

png("plots/revised_cisLASSO.png", width = 13.5, height = 9, units = "in", res = 300)
prot_cisresid_models$p_left + prot_cisresid_models$p_right
dev.off()



#---------------------------------
### plot calibration error
#---------------------------------

ece= ggplot(prot_resid_models$df_wide,
            aes(x = `ECE_genetically adjusted proteins`,
                y = `ECE_unadjusted proteins`)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific) +
  coord_cartesian(xlim = c(0, 0.0125), ylim = c(0, 0.012)) +
  
  geom_text_repel(
    aes(
      label = ifelse((abs(`ECE_genetically adjusted proteins`-`ECE_unadjusted proteins`))>0.001,
        Endpoint,
        NA)
    ),
    size = 3,
    hjust = 1,
    max.overlaps = 10,
    box.padding = 0.4,
    point.padding = 0.3,
    
  ) + 
  theme_classic() +
  labs(
    x = "ECE (genetically adjusted proteins)",
    y = "ECE (unadjusted proteins)"
  ) +
  theme(axis.title.x = element_text(face = "bold",size = 12),
        axis.title.y = element_text(face = "bold",size = 12),
        axis.text.x = element_text(size = 11),   # <-- tick labels
        axis.text.y = element_text(size = 11),
        plot.tag = element_text(face = "bold", size = 14))
  

all= readRDS("data/revised_LASSO_AUC_all_models_extracted.rds")

cal= rbind(subset(all[[2]]$predictions, Disease=="KNEE_ARTHROSIS_DATE"),
           subset(all[[1]]$predictions, Disease=="KNEE_ARTHROSIS_DATE"))

c= calibration_plot(data = cal,
                 obs = "y", pred = "pred",
                 group = "Model",
                 title = "Knee-Osteoarthritis")$calibration_plot + theme(plot.tag = element_text(face = "bold", size = 14),
                                                                            legend.position = c(0.75, 0.2),
                                                                         axis.title.x = element_text(face = "bold",size = 12),
                                                                         axis.title.y = element_text(face = "bold",size = 12),
                                                                         axis.text.x = element_text(size = 11),   # <-- tick labels
                                                                         axis.text.y = element_text(size = 11))
pdf("plots/AUC_calibration.pdf", width = 15.5, height= 7)
ece + c + plot_annotation(tag_levels = 'A')
dev.off()

png("plots/AUC_calibration.png", width = 15.5, height= 7, units = "in", res = 300)
ece + c + plot_annotation(tag_levels = 'A')
dev.off()



### C-index ECE 5 
ece5= readRDS("data/ECE5_cindex_plot.rds")

pdf("plots/calibrations_comparison.pdf", width = 17, height= 8)
ece + ece5 +  plot_annotation(tag_levels = 'A')
dev.off()


png("plots/calibrations_comparison.png", width = 17, height= 8, units = "in", res = 300)
ece + ece5 +  plot_annotation(tag_levels = 'A')
dev.off()
