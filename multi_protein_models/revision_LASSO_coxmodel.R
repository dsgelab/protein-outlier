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
library(survcomp, lib = "~/R_library/")
library(ggrepel)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")


## Run stability selection when selecting features with LASSO - cox-model and c-index extraction  - And calibration 
lasso_cox_with_stability <- function(X, Y, time_mat, status_mat, seed = 1234, testsize = 0.3,
                                 cutoff = 0.6, PFER = 3, B = 100,
                                 nfolds_cv = 5, nlambda = 200,
                                 ece_bins = 15,    # number of bins for ECE/ICI
                                 years= 5,  ## followup years to calculate calibration
                                 ncores=15) {
  diseases <- colnames(time_mat)
  
  # ---- output containers ----
  results <- data.frame(
    Disease   = diseases,
    C_index       = NA_real_,
    C_Lower = NA_real_,
    C_Upper = NA_real_,
    n_features  = NA_integer_,
    features    = NA_character_,
    # calibration / proper scoring
    ECE5        = NA_real_,
    n_eval5   = NA_integer_,
    n_events5 = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  pred_list <- vector("list", length(diseases))
  
  # ---- parallel backend ----
  #cl <- makeCluster(parallelly::availableCores(omit = 2))
  #registerDoParallel(cl)
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  # always stop cluster even if function errors
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  # avoid nested foreach inside workers
  parallel::clusterEvalQ(cl, {
    library(foreach)
    foreach::registerDoSEQ()
    NULL
  })
  
  out <- foreach(i = seq_along(diseases),
                 .packages = c("caret","glmnet","pROC","stats","survival"),
                 .errorhandling = "pass") %dopar% {
                   
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
                   library(survcomp)

                   set.seed(seed)
                   
                   # ---- outcome ----
                   # ---- outcome (keep IDs) ----
                   y0 <- Y[, i]
                   ids <- names(y0)          
                   ids <- names(y0)[!is.na(y0)]    # keep non-missing outcome
                   
                   # align X, time, status to the SAME ids and the SAME disease column i
                   X_filt <- X[ids, , drop = FALSE]
                   tt <- time_mat[ids, i]
                   dd <- status_mat[ids, i]
                   
                   # drop any remaining missingness jointly 
                   keep <- is.finite(tt) & is.finite(dd)
                   ids2 <- ids[keep]
                   
                   y  <- as.integer(y0[ids2])
                   X_filt <- X_filt[ids2, , drop = FALSE]
                   tt <- as.numeric(tt[keep])
                   dd <- as.integer(dd[keep])
                  

                   # ---- train/test split ----
                   train_idx <- caret::createDataPartition(y, p = 1 - testsize, list = FALSE)
                   
                   y_train <- y[train_idx] ### endpoint splitting
                   y_test  <- y[-train_idx]
                   
                   X_train <- X_filt[train_idx, , drop = FALSE] ### proteins splitting
                   X_test  <- X_filt[-train_idx, , drop = FALSE]
                   
                   ### TIME AND STATUS MATRICES splitting
                   t_train <- tt[train_idx] 
                   t_test  <- tt[-train_idx]
                   
                   d_train <- dd[train_idx] 
                   d_test  <- dd[-train_idx]
                   
                   
                   # ---- stability selection on TRAIN only USING LOGISTIC REGRESSION ----
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
                     piece <- list(
                       C_index = NA_real_,
                       C_Lower= NA_real_,
                       C_Upper= NA_real_,
                       final_vars = character(0),
                       ECE5 = NA_real_,
                       n_eval5 = NA_integer_,
                       n_events5 = NA_integer_,
                       p5_test = rep(NA_real_, length(t_test)),
                       known5  = NA_integer_,  # full length = length(d_test)
                       y5      = NA_integer_,      # full length = length(d_test)
                       y5k     = NA_integer_,     # subset
                       p5k     = NA_integer_,     # subset
                       t_test = NA_integer_,
                       d_test = NA_integer_,
                       stable_vars = stable_vars,
                       selprob = stab_fit$max
                     )
                     return(list(i = i, piece = piece))
                   }

                   
                   # ---- SPECIAL CASE: only 1 stable feature 
                   if (length(stable_idx) == 1L) {
                     
                     x1_train <- as.numeric(X_train[, stable_idx])
                     x1_test  <- as.numeric(X_test[,  stable_idx])
                     
                     # estimate weight for this single feature on TRAIN (unpenalized)
                     #fit1 <-  survival::coxph(survival::Surv(t_train, d_train) ~ x1_train)
                     #w_final <- as.numeric(stats::coef(fit1)[1]) # slope for the feature - note that cox model has no intercept so it is [1]
                     
                     fit1 <- glm(y_train ~ x1_train, family = "binomial")
                     w_final <- as.numeric(coef(fit1)[2])          # slope for the feature
                     
                     # build score (single-feature score)
                     score_train <- x1_train * w_final
                     score_test  <- x1_test  * w_final
                     
                     final_vars = colnames(X_train)[stable_idx]
                     
                     data= data.frame(score= score_train,
                                      time = t_train,
                                      status= d_train)
                     
                   } else {
                     
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
                     sel2 <- which(w2 != 0) # take non zero features
                     
                     if (length(sel2) == 0L) {
                       # If lambda.min shrinks all to 0, try lambda.1se as a conservative alternative
                       coefs2 <- coef(cv_fit2, s = "lambda.1se")
                       w2 <- as.numeric(coefs2[-1, , drop = FALSE])
                       sel2 <- which(w2 != 0)
                     }
                     
                     # if still none selected, return NA
                     if (length(sel2) == 0L) {
                       piece <- list(
                         C_index = NA_real_,
                         C_Lower= NA_real_,
                         C_Upper= NA_real_,
                         final_vars = character(0),
                         ECE5 = NA_real_,
                         n_eval5 = NA_integer_,
                         n_events5 = NA_integer_,
                         p5_test = rep(NA_real_, length(t_test)),
                         known5  = NA_integer_,  # full length = length(d_test)
                         y5      = NA_integer_,      # full length = length(d_test)
                         y5k     = NA_integer_,     # subset
                         p5k     = NA_integer_,     # subset
                         t_test = NA_integer_,
                         d_test = NA_integer_,
                         stable_vars = stable_vars,
                         selprob = stab_fit$max
                       )
                       return(list(i = i, piece = piece))
                     }
                     
                     # map sel2 (within stable set) back to original column indices
                     cols_final <- stable_idx[sel2]
                     w_final <- w2[sel2]
                     final_vars <- colnames(X_train)[cols_final]
                     
                     #---- build linear score 
                     score_train <- as.numeric(X_train[, cols_final, drop = FALSE] %*% w_final)
                     score_test  <- as.numeric(X_test[, cols_final, drop = FALSE] %*% w_final)
                     
                     data= data.frame(score= score_train,
                                      time = t_train,
                                      status= d_train)
                     
                   }
                     
                     ### fit score on training set
                     model <- survival::coxph(survival::Surv(time, status) ~ score, data= data)
                     
                     ## predict on TEST set 
                     pred <- predict(model,
                                     newdata = data.frame(score = score_test),
                                     type = "lp")
                     
                     # ---- C-index on TEST ----
                     
                     # d_test should be 0/1 with 1 = event
                     ci <- survcomp::concordance.index(
                       x = pred,
                       surv.time = t_test,
                       surv.event = d_test,
                       method = "noether"   # common variance estimator
                     )
                     
                     C_index <- as.numeric(ci$c.index)
                     C_Lower <- as.numeric(ci$lower)
                     C_Upper <- as.numeric(ci$upper)
                     
                     # ---- Predicted 5-year risk on TEST ----
                     # baseline cumulative hazard from TRAIN model
                     
                     # horizon in SAME units as t_test / bh$time
                     horizon <- years * 365.25
                     
                     bh <- survival::basehaz(model, centered = FALSE)
                     H0t <- stats::approx(bh$time, bh$hazard, xout = horizon, rule = 2)$y
                     S0t <- exp(-H0t)
                     
                     # S(e|x) = S0(5) ^ exp(lp)
                     p5_test <- 1 - (S0t ^ exp(pred))
                     
                     # ---- 5-year calibration subset (exclude censored before 5y) ----
                     known5 <- (d_test == 1L & t_test <= horizon) | (d_test == 0L & t_test >= horizon)
                     y5 <- ifelse(d_test == 1L & t_test <= horizon, 1L, 0L)
                     
                     y5k <- y5[known5]
                     p5k <- p5_test[known5]
                     
                     n_eval5 <- length(y5k)
                     n_events5 <- sum(y5k == 1L)
                     
                     ECE5 <- if (length(unique(y5k)) == 2L) {
                       ece_binary(y5k, p5k, n_bins = ece_bins)
                     } else {
                       NA_real_
                     }
                     
                     piece <- list(
                       C_index = as.numeric(C_index),
                       C_Lower = C_Lower,
                       C_Upper = C_Upper,
                       final_vars = final_vars,
                       ECE5 = ECE5,
                       n_eval5 = n_eval5,
                       n_events5 = n_events5,
                       p5_test = p5_test,
                       known5  = known5,  # full length = length(d_test)
                       y5      = y5,      # full length = length(d_test)
                       y5k     = y5k,     # subset
                       p5k     = p5k,     # subset
                       t_test = t_test,
                       d_test = d_test,
                       stable_vars = stable_vars,
                       selprob = stab_fit$max
                     )
                     
                     list(i = i, piece = piece)
                   
                 }
                   
                   parallel::stopCluster(cl)
                   
                   
                   # ---- collect results ----
                   for (obj in out) {
                     if (is.null(obj) || inherits(obj, "error") || is.null(obj$piece)) next
                     
                     i <- obj$i
                     piece <- obj$piece
                     
                     results$C_index[i] <- piece$C_index
                     results$C_Lower[i] <- piece$C_Lower
                     results$C_Upper[i] <- piece$C_Upper
                     results$ECE5[i]   <- piece$ECE5
                     results$n_eval5[i]   <- piece$n_eval5
                     results$n_events5[i] <- piece$n_events5
                     
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

time_df   = data.frame(ID = intervene_fi$ID)
status_df = data.frame(ID = intervene_fi$ID)
all_endpoints= data.frame(ID = intervene_fi$ID)

for (dis in colnames(intervene_fi %>% select(contains("_DATE")))) {
  df= intervene_fi %>% select(ID, END_OF_FOLLOWUP, dis, baseline) 
  df_merge= df %>% mutate(fol = difftime(.data[[dis]], baseline, units = "weeks")/52.25) ###note that I converted in years!!
  
  ## excluding all prevalent cases or incident cases recorded within the first 6 months of follow-up
  df_filt= df_merge %>% filter(fol>0.5 | is.na(fol))
  
  ### binary outcome, 1= incidence, 0 = right censored 
  ## MATRIX FOR EXTRACTING FEATURES WITH LOGIF
  df_bin= df_filt %>% mutate(!!as.name(dis) := ifelse(is.na(!!as.name(dis)),0,1)) %>% select(ID, all_of(dis))
  
  df_filt = df_filt %>%
    mutate(time_to_event = as.numeric(ifelse(is.na(.data[[dis]]),
                                             (difftime(END_OF_FOLLOWUP, baseline, units="days")),  ### end of followup - baseline
                                             (difftime(.data[[dis]], baseline, units = "days")))), ### incidence - baseline 
           status = ifelse(is.na(fol), 0, 1)  ### 1= incidence , 0= health control
    ) %>%
    select(ID, time_to_event, status)
  
  all_endpoints = full_join(all_endpoints, df_bin, by = "ID")
  
  ### MATRICES FOR COX MODELS 
  time_df = time_df %>%
    full_join(
      df_filt %>% dplyr::select(ID, time_to_event),
      by = "ID"
    ) %>%
    rename(!!dis := time_to_event)
  
  status_df = status_df %>%
    full_join(
      df_filt %>% dplyr::select(ID, status),
      by = "ID"
    ) %>%
    rename(!!dis := status)
}

# set rownames = ID, drop ID column, convert to matrices
rownames(time_df) = time_df$ID
time_df$ID= NULL
m_time= as.matrix(time_df)

rownames(status_df)= status_df$ID
status_df$ID= NULL
m_status= as.matrix(status_df)

rownames(all_endpoints)= all_endpoints$ID
all_endpoints$ID= NULL
m_allendpoints= as.matrix(all_endpoints)

### scale, and create matrix
load_scale_filter_matrix <- function(file, match_ids, id_col = 1) {
  dt <- fread(file)
  setDT(dt)[, (names(dt)[-id_col]) := lapply(.SD, scale), .SDcols = names(dt)[-id_col]]
  
  mat <- as.matrix(dt, rownames = TRUE)
  mat <- mat[rownames(mat) %in% match_ids, ]
  return(mat)
}

#------------------------------------------------ LOAD ALL THE PREDICTORS ---------------------------------------
m_allpheno= load_scale_filter_matrix("data/idefix_olinksoma/pheno_residuals.tsv", rownames(m_time))

m_resid= load_scale_filter_matrix("data/idefix_olinksoma/residuals_output_nocov.tsv", rownames(m_time))
# Drop column after loading
m_resid= m_resid[, colnames(m_resid) != "AggAbsResid"]


m_cisresid= load_scale_filter_matrix("data/cis_residuals_output_nocov.tsv", rownames(m_time))
m_cisresid= m_cisresid[, colnames(m_cisresid) != "AggAbsResid"]



##------------------------------ RUN LASSO ----------------------------------------------------------------------------

#################### run LASSO   

### proteins
pred_prot= lasso_cox_with_stability(m_allpheno,m_allendpoints,m_time,m_status)

### residuals
pred_resid= lasso_cox_with_stability(m_resid,m_allendpoints,m_time,m_status)

######### extract dataframes
extract_results <- function(model_result, model_name) {
  df <- as.data.frame(model_result[[1]])
  colnames(df) <- gsub("^results\\.", "", colnames(df))
  df$Model <- model_name
  
  preds <- model_result[[2]]
  pred_long <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    pi <- preds[[i]]
    if (is.null(pi) || is.null(pi$d_test) || is.null(pi$p5_test)) next
    
    n <- length(pi$d_test)
    
    # counts (test-set)
    df$n_total[i]    <- n
    df$n_cases[i]    <- sum(pi$d_test == 1L, na.rm = TRUE)
    df$n_controls[i] <- sum(pi$d_test == 0L, na.rm = TRUE)
    
    pred_long[[i]] <- data.frame(
      Disease = rep(df$Disease[i], n),
      Model   = rep(model_name, n),
      t_test  = pi$t_test,
      d_test  = pi$d_test,
      p5_test = pi$p5_test,
      known5  = if (!is.null(pi$known5)) pi$known5 else rep(NA, n),
      y5      = if (!is.null(pi$y5))     pi$y5     else rep(NA, n)
    )
  }
  
  list(
    summary = df,
    predictions = dplyr::bind_rows(pred_long)
  )
}


all_results <- list(
  proteins        = extract_results(pred_prot, "unadjusted proteins"),
  resid           = extract_results(pred_resid, "genetically adjusted proteins")
)

saveRDS(all_results, file = "data/revision_LASSO_coxmodel_predictions_all.rds")


lasso= readRDS("data/revision_LASSO_coxmodel_predictions_all.rds")

# 1) combined summary table across models
summary_df <- bind_rows(lapply(lasso, `[[`, "summary"))

# 2) combined individual-level prediction table across models
pred_df <- bind_rows(lapply(lasso, `[[`, "predictions"))

feature_table <- summary_df %>%
  select(Disease, Model, n_features) %>%
  pivot_wider(
    names_from  = Model,
    values_from = n_features
  )

### All predictions
summary_df$Disease= gsub("_DATE","",summary_df$Disease)
summary_df= merge(summary_df,definitions,
               by.x="Disease", by.y="FinnGen endpoint")
summary_df$features= gsub("StandResid_","",summary_df$features)

write.table(summary_df, "data/revision_LASSO_Cindex_summarydf.tsv",quote=F, sep="\t",row.names = F)


allpred= fread("data/revision_LASSO_Cindex_summarydf.tsv")
### standard error
allpred= allpred %>% mutate(SE= (C_Upper-C_Lower)/3.92)


# ---- choose models once ----

### data-frame , metrics and plot between 2 Models 
compare_c_models <- function(allpred, models_keep, compare = models_keep) {
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
      values_from = c(C_index, C_Lower, C_Upper, SE, ECE5, n_features, features),
      names_glue = "{.value}_{Model}"
    )
  
  # dynamic cols
  C1 <- paste0("C_index_", m1); C2 <- paste0("C_index_", m2)
  SE1  <- paste0("SE_",  m1); SE2  <- paste0("SE_",  m2)
  
  # ---- z-test (one-tailed: m2 > m1) ----
  df_wide <- df_wide %>%
    mutate(
      zscore = (.data[[C2]] - .data[[C1]]) / sqrt(.data[[SE2]]^2 + .data[[SE1]]^2),
      pval_z = 1 - pnorm(zscore),
      sig_zscore = if_else(pval_z < 0.05, "yes", "no"),
      diff = .data[[C2]] - .data[[C1]],
      perc_gain = (diff / .data[[C1]]) * 100
    )
  
  # ---- join back for plotting ----
  allpred_out <- allpred_f %>%
    left_join(
      df_wide %>% select(Disease, Endpoint, sig_zscore),
      by = c("Disease", "Endpoint")
    ) %>%
    mutate(
      C_Lower = pmax(C_Lower, 0.5),
      C_index       = pmax(C_index, 0.5),
      alpha_ztest = if_else(sig_zscore == "yes", 1, 0.5, missing = 0.5)
    )
  
  
  endpoint_levels <- allpred_out %>%
    group_by(Endpoint) %>%
    summarise(Max_C = max(C_index, na.rm = TRUE), .groups = "drop") %>%
    arrange(Max_C) %>%              # low -> high
    pull(Endpoint)
  
  allpred_out <- allpred_out %>%
    group_by(Endpoint) %>%
    mutate(Max_C = max(C_index, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      Endpoint = reorder(Endpoint, Max_C),
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
    y = Endpoint, x = C_index,
    xmin = C_Lower, xmax = C_Upper,
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
    labs(x = "C-index (95% CI)", y = "Disease") +
    annotate(
      "text", y = 1.5, x = 0.7,
      label = expression("* One-tailed p < 0.05 ("*Delta*"C-index z-score)"),
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

prot_resid_models= compare_c_models(allpred, c("unadjusted proteins","genetically adjusted proteins"))

ece= ggplot(prot_resid_models$df_wide,
            aes(x = `ECE5_genetically adjusted proteins`,
                y = `ECE5_unadjusted proteins`)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific) +
  coord_cartesian(xlim = c(0, 0.0125), ylim = c(0, 0.012)) +
  
  geom_text_repel(
    aes(
      label = ifelse((abs(`ECE5_genetically adjusted proteins`-`ECE5_unadjusted proteins`))>0.0005,
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
    x = "5-year risk ECE (genetically adjusted proteins)",
    y = "5-year risk ECE (unadjusted proteins)"
  ) +
  theme(axis.title.x = element_text(face = "bold",size = 12),
        axis.title.y = element_text(face = "bold",size = 12),
        axis.text.x = element_text(size = 11),   # <-- tick labels
        axis.text.y = element_text(size = 11),
        plot.tag = element_text(face = "bold", size = 14))

saveRDS(ece, "data/ECE5_cindex_plot.rds")


############ 
pdf("plots/revised_coxLASSO.pdf", width = 11, height = 10)
prot_resid_models$p_left
dev.off()

png("plots/revised_coxLASSO.png", width = 11, height = 10, units = "in", res= 300)
prot_resid_models$p_left
dev.off()


