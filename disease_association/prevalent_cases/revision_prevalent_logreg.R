library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggrepel)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

## function to scale dataframe
scale_df <- function(df, idcol) {
  cols <- setdiff(names(df), idcol)
  df[, (cols) := lapply(.SD, scale), .SDcols = cols]
  df
}

### FUNCTION FOR  PREVALENT LOGISTIC REGRESSION FOR FILTERED PAIRS #################
prevalent_logif_filtpairs <- function(prot_df, idcol,
                                  endpoint_date_suffix = "_DATE",
                                  sep = "__") {
  
  stopifnot(idcol %in% names(prot_df))
  
  date_cols <- names(intervene_fi %>% select(contains(endpoint_date_suffix)))
  
  res_all <- list()
  
  for (score_col in setdiff(names(prot_df), idcol)) {
    
    parts <- str_split_fixed(score_col, fixed(sep), 2)
    prot  <- parts[1]
    endp  <- parts[2]
    
    dis_date <- paste0(endp, endpoint_date_suffix)
    
    df_base <- intervene_fi %>%
      select(ID, baseline, age, sex, all_of(dis_date))
    
    p <- prot_df %>%
      select(ID = all_of(idcol), score = all_of(score_col))
    
    df_merge <- merge(df_base, p, by = "ID")
    
    ### binary outcome, 1 = prevalent at baseline, 0 = not prevalent
    df_bin <- df_merge %>%
      mutate(y = ifelse(!is.na(.data[[dis_date]]) &
                          .data[[dis_date]] <= baseline, 1, 0))
    
    fit <- glm(y ~ score + age + sex, data = df_bin, family = "binomial")
    
    out = as.data.frame(summary(fit)$coefficients["score", , drop = FALSE])
    out$protein  <- prot
    out$endpoint <- gsub("_DATE","",dis_date)
    out$n        <- nrow(df_bin)
    out$cases    <- sum(df_bin$y == 1, na.rm = TRUE)
    
    res_all[[length(res_all) + 1]] = out
  }
  
  bind_rows(res_all)
}

## FUNCTION FOR PREVALENT LOGISTIC REGRESSION FOR ALL PROTEINS AND RESIDUALS NON FILTERED
prevalent_logif= function(prot_df, idcol) {
  allendpoints = data.frame()
  for (dis in colnames(intervene_fi %>% select(contains("_DATE")))) {
    df= intervene_fi %>% select(ID, END_OF_FOLLOWUP, dis, baseline, age, sex) 
    allresults= data.frame()
    for (prot in colnames(prot_df)[-1]) {
      p= prot_df %>% select(ID= idcol, prot)
      df_merge= merge(df, p, by="ID")
      
      ### binary outcome, 1 = prevalent at baseline, 0 = not prevalent
      df_bin= df_merge %>% mutate(!!as.name(dis) := ifelse(!is.na(!!as.name(dis)) & as.Date(!!as.name(dis)) <= as.Date(baseline), 1, 0))
      
      ### Logistic regression
      formula= as.formula(paste0(dis,"~",prot," + age + sex"))
      glm_model = glm(formula, data = df_bin, family = "binomial")
      
      glmdf= as.data.frame(summary(glm_model)$coefficients[2, , drop = FALSE])
      glmdf$protein= prot
      
      allresults = rbind(glmdf, allresults)
      
    }
    allresults$endpoint= dis
    allendpoints= rbind(allresults, allendpoints)
  }
  return(allendpoints)
}


#------- LOAD DATASETS ----------------------------------

residuals_filt_cis= fread("data/filtered_cis_pPGS/cis_residuals_output_nocov.tsv") %>% select(!(AggAbsResid)) 
###scale 
residuals_filt_cis = scale_df(residuals_filt_cis, idcol = "IndID")

residuals_filt_trans= fread("data/filtered_trans_pPGS/trans_residuals_output_nocov.tsv") %>% select(!(AggAbsResid))
###scale 
residuals_filt_trans = scale_df(residuals_filt_trans, idcol = "IndID")

pgs_filt_cis= fread("data/filtered_cis_pPGS/cis_genoresiduals.tsv") 
###scale 
pgs_filt_cis = scale_df(pgs_filt_cis, idcol = "eid")

pgs_filt_trans= fread("data/filtered_trans_pPGS/trans_genoresiduals.tsv") 
###scale 
pgs_filt_trans = scale_df(pgs_filt_trans,  idcol = "eid")

residuals_filt_gw= fread("data/filtered_gw_pPGS/gw_residuals_output_nocov.tsv") %>% select(!(AggAbsResid)) 
###scale 
residuals_filt_gw = scale_df(residuals_filt_gw,  idcol = "IndID")

pgs_filt_gw= fread("data/filtered_gw_pPGS/gw_genoresiduals.tsv") 
###scale 
pgs_filt_gw = scale_df(pgs_filt_gw,  idcol = "eid")


allproteins=fread("data/idefix_olinksoma/pheno_residuals.tsv")
###scale proteins
allproteins= scale_df(allproteins, idcol="eid")

##### residuals with no adjustment for age and sex
resid= fread("data/idefix_olinksoma/residuals_output_nocov.tsv")
resid$AggAbsResid=NULL
##scale 
resid= scale_df(resid, idcol="IndID")

pgs= fread("data/idefix_olinksoma/geno_residuals.tsv")
###scale pgs
pgs= scale_df(pgs, idcol="eid")

### file with first occurrences of INTERVENE phenotypes 
intervene_fi= arrow::read_parquet("/scratch/project_2007428/projects/prj_010_phenotype_file/share/ukb78537_INTERVENE_phenotype-2024-08-01.parquet")
intervene_fi= intervene_fi %>% select(ID,SEX,END_OF_FOLLOWUP, contains("_DATE"))
### filter for our individuals
intervene_fi = intervene_fi %>% filter(ID %in% residuals_filt_cis$IndID)

### filter phenotypes for the intervene. note that from the 39 endpoints, BMI and Covid-19 severity are not present as first occurrences, so the final number is 37 
definitions= fread("data/Intervene_flagship_endpoint_collection_Definitions.tsv") 
definitions = definitions %>% select(Endpoint,`FinnGen endpoint`) 

## load age and sex covariates, first assessment visit date and death
phenos=fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv")

selphenos= phenos %>% select(ID=eid, age='21022-0.0',
                             sex='31-0.0', baseline='53-0.0', death='40000-0.0')
selphenos = selphenos %>% filter(ID %in% intervene_fi$ID)

#### The column "END_OF_FOLLOWUP" is the age of the first record of a disease diagnosis (for individuals with the diseases), age at death for a cause other than the disease, age at last record available in the registries or electronic health records or age 80, whatever happened first
intervene_fi= merge(intervene_fi, selphenos, by="ID")


#### -------- RUN PREVALENT CASES ON ORIGINAL PAIRS ----------------

prev_proteins= prevalent_logif(prot_df= allproteins, idcol = "eid")
prev_proteins$Model= "protein"

prev_residuals= prevalent_logif(prot_df= resid, idcol = "IndID")
prev_residuals$Model= "residuals"

prev_pgs= prevalent_logif(prot_df= pgs, idcol = "eid")
prev_pgs$Model= "pPGS"


### merge all 
all= bind_rows(prev_proteins, prev_residuals, prev_pgs)

all = all %>% rename(zscore = `z value`,
                     pval= "Pr(>|z|)",
                     se= "Std. Error")

all$protein= gsub("StandResid_","",all$protein)

write.table(all, "data/revision_prevalent_logisticmodels_protdisease.tsv", quote=F, sep="\t", row.names = F)






