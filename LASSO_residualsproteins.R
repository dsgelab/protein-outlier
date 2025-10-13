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
library(PRROC, lib.loc="~/R_library/")
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")


#### fit the lasso and extract AUC 
lasso= function(X, Y, seed = 1234, testsize=0.3) {
  
  results <- data.frame(
    Disease = colnames(Y),
    PR_AUC= NA,
    PR_Lower=NA,
    PR_Upper=NA,
    AUC = NA,
    AUC_Lower = NA,
    AUC_Upper = NA
  )
  
  pred_list <- list()
  
  # register parallel backend
  cl = makeCluster(parallelly::availableCores(omit = 1))
  registerDoParallel(cl)
  
  out = foreach(i = 1:ncol(Y), .packages = c("caret","glmnet","pROC"),
                .errorhandling = "remove") %dopar% {
    .libPaths("~/R_library/") 
    library("PRROC")             
                  
    
    set.seed(seed)  
    
    y = Y[, i] ## extract the i-th outcome
    y= y[!is.na(y)] ## Remove NAs (which are excluded individuals, prevalent cases)
    if (length(unique(y)) <= 1) return(NULL)
    
    X_filt= X[names(y),]
    
    #### splitting data in training and test
    train_idx= createDataPartition(y, p = 1 - testsize, list = FALSE)    
    y_train= y[train_idx] ### endpoint splitting
    y_test= y[-train_idx]
    
    X_filt_train= X_filt[train_idx, ] ### proteins splitting
    X_filt_test= X_filt[-train_idx, ]
    
    cv_fit = cv.glmnet(X_filt_train, y_train,
                       family = "binomial",
                       alpha = 1,
                       #lambda.min.ratio = 1e-4,
                       nfolds = 10,
                       nlambda= 200,
                       type.measure="auc") #Fit a LASSO logistic regression model using cross-validation
    
    coefs= coef(cv_fit, s = "lambda.min") #Store coefficients at best lambda 
    
    w=as.numeric(coefs[-1, , drop = FALSE]) # Drop intercept, get weights as numeric vector
    sel_prot= which(w != 0) # Identify non-zero (selected) features

    if (length(sel_prot) == 0) {
      # If no features are selected, set to NA
      results$AUC[i] = NA
      results$AUC_Lower[i] = NA
      results$AUC_Upper[i] = NA
      next  # Skip to next outcome
    }
    # Compute score as linear combination of selected features (weighted sum)
    scores_prot_train= X_filt_train[, sel_prot, drop = FALSE] %*% w[sel_prot]
    ### scale
    scores_prot_train= scale(scores_prot_train)
    
    # Compute score as linear combination of selected features (weighted sum)
    scores_prot_test= X_filt_test[, sel_prot, drop = FALSE] %*% w[sel_prot]
    ### scale
    scores_prot_test= scale(scores_prot_test)
    
    train_df= data.frame(
      y = y_train,
      score = scores_prot_train
    )
    
    ### train the model
    model= glm(y ~ score, data=train_df,family = "binomial")
    
    test_df= data.frame(
      y= y_test,
      score = scores_prot_test
    )
    
    ### predict on test data
    pred= predict(model, newdata= test_df, type = "response")
    
    ### extract also PR - AUC
    # --- PR AUC + CI with PRROC ---
    yt <- as.integer(y_test) # ensure numeric 0/1
    fg <- pred[yt == 1L]    # scores for positives
    bg <- pred[yt == 0L]    # scores for negatives
    
    # point estimate
    pr_obj  <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    
    prauc  <- as.numeric(pr_obj$auc.integral)
    
    # bootstrap CI (stratified)
    B <- 200L                                   
    set.seed(seed + i)                          # reproducible per-outcome
    if (length(fg) >= 2L && length(bg) >= 2L) {
      npos <- length(fg); nneg <- length(bg)
      aucs <- replicate(B, {
        sfg <- sample(fg,  npos, replace = TRUE)
        sbg <- sample(bg,  nneg, replace = TRUE)
        PRROC::pr.curve(scores.class0 = sfg, scores.class1 = sbg, curve = T)$auc.integral
      })
      prlower <- as.numeric(quantile(aucs, 0.025, na.rm = TRUE))
      prupper <- as.numeric(quantile(aucs, 0.975, na.rm = TRUE))
    } else {
      prlower <- NA_real_; prupper <- NA_real_
    }
    
    res = c(
      PR_AUC   = prauc,
      PR_Lower = prlower,
      PR_Upper = prupper,
      AUC       = auc(y_test, pred),
      AUC_Lower = ci.auc(y_test, pred)[1],
      AUC_Upper = ci.auc(y_test, pred)[3]
    )
    
    list(
      i   = i,
      res = res,
      pred = list(
        disease = colnames(Y)[i],
        y_test = y_test,
        pred = pred,
        pr_obj = pr_obj
      )
    )
    
  }
  
  stopCluster(cl)
  
  # collect results back into same format
  for (piece in out) {
    if (is.null(piece) || is.null(piece$res)) next
    i <- piece$i
    results$PR_AUC[i] <- piece$res["PR_AUC"]
    results$PR_Lower[i] <- piece$res["PR_Lower"]
    results$PR_Upper[i] <- piece$res["PR_Upper"]
    results$AUC[i]       <- piece$res["AUC"]
    results$AUC_Lower[i] <- piece$res["AUC_Lower"]
    results$AUC_Upper[i] <- piece$res["AUC_Upper"]
    pred_list[[i]]       <- piece$pred
  }
  
  return(list(results = results, predictions = pred_list))
}


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
m_allprs= load_scale_filter_matrix("data/idefix_olinksoma/geno_residuals.tsv", rownames(m_allendpoints))
m_residnocov= load_scale_filter_matrix("data/idefix_olinksoma/residuals_output_nocov.tsv", rownames(m_allendpoints))
# Drop column after loading
m_residnocov= m_residnocov[, colnames(m_residnocov) != "AggAbsResid"]

#### proteins + prs
m_protprs= merge(m_allpheno,m_allprs,by="row.names", suffixes= c("prot_","prs_"))
rownames(m_protprs)= m_protprs$Row.names
m_protprs$Row.names=NULL
m_protprs= as.matrix(m_protprs)    ### matrix PRS +proteins 

m_cisprs= load_scale_filter_matrix("data/cis_genoresiduals.tsv", rownames(m_allendpoints))

####### combined protein + cis-prs 
m_protcisprs= merge(m_allpheno,m_cisprs,by="row.names", suffixes= c("prot_","cisprs_"))
rownames(m_protcisprs)= m_protcisprs$Row.names
m_protcisprs$Row.names=NULL
m_protcisprs= as.matrix(m_protcisprs)    ### matrix PRS +proteins 

##------------------------------ RUN LASSO ----------------------------------------------------------------------------

#################### run LASSO NO AGE AND SEX   

### proteins
pred_prot= lasso(m_allpheno,m_allendpoints)

### residuals
pred_resid_nocov= lasso(m_residnocov, m_allendpoints)

### PRS
pred_prs= lasso(m_allprs, m_allendpoints)

#### Combined proteins + PRS
prot_prs= lasso(m_protprs,m_allendpoints)

##### cis-PRS
pred_cisprs= lasso(m_cisprs,m_allendpoints)

### combined proteins + cis-PRS
pred_protcisprs= lasso(m_protcisprs,m_allendpoints)



######### extract dataframes
extract_auc_results <- function(model_result, model_name) {
  df <- as.data.frame(model_result[[1]])
  colnames(df) <- gsub("^results\\.", "", colnames(df))
  df$Model <- model_name
  return(df)
}

prot_auc = extract_auc_results(pred_prot, "unadjusted proteins")
resid_auc = extract_auc_results(pred_resid_nocov, "genetically adjusted proteins")

prs_auc = extract_auc_results(pred_prs, "PGS")

protprs_auc = extract_auc_results(prot_prs, "unadjusted proteins-PGS combined")
cisprs_auc= extract_auc_results(pred_cisprs, "cis-PGS")
protcisprs_auc= extract_auc_results(pred_protcisprs, "unadjusted proteins-cis-PGS combined")


### All predictions
allpred= rbind(prot_auc, resid_auc, prs_auc,protprs_auc, cisprs_auc, protcisprs_auc)
allpred$Disease= gsub("_DATE","",allpred$Disease)
allpred= merge(allpred,definitions,
               by.x="Disease", by.y="FinnGen endpoint")
write.table(allpred, "data/LASSO_polyprotein_all.tsv",quote=F, sep="\t",row.names = F)


allpred= rbind(prot_auc,resid_auc)
allpred$Disease= gsub("_DATE","",allpred$Disease)
allpred= merge(allpred,definitions,
               by.x="Disease", by.y="FinnGen endpoint")
#write.table(allpred, "data/LASSO_polyprotein_protresid.tsv",quote=F, sep="\t",row.names = F)



#allpred= fread("data/LASSO_polyprotein_protresid.tsv")
### standard error
allpred= allpred %>% mutate(SE= (AUC_Upper-AUC_Lower)/3.92)

df_wide= allpred %>% pivot_wider(
  id_cols = c(Disease, Endpoint),
  names_from = Model,
  values_from = c(AUC, AUC_Lower, AUC_Upper, SE),
  names_glue = "{.value}_{Model}") 

### z-score of the difference and one-tailed p value
df_wide= df_wide %>% mutate(zscore= (`AUC_genetically adjusted proteins`-`AUC_unadjusted proteins`)/sqrt(`SE_genetically adjusted proteins`^2+`SE_unadjusted proteins`^2),
                            pval_z= (1 - pnorm(abs(zscore)))
                            ) 

