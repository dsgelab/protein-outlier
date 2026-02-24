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

residuals_filt_cis= fread("data/FG_filtered_cis_pPGS/cis_residuals_output_nocov.tsv") %>% select(!(AggAbsResid)) 
###scale 
residuals_filt_cis = scale_df(residuals_filt_cis, idcol = "IndID")

residuals_filt_trans= fread("data/FG_filtered_trans_pPGS/trans_residuals_output_nocov.tsv") %>% select(!(AggAbsResid))
###scale 
residuals_filt_trans = scale_df(residuals_filt_trans, idcol = "IndID")

pgs_filt_cis= fread("data/FG_filtered_cis_pPGS/cis_genoresiduals.tsv") 
###scale 
pgs_filt_cis = scale_df(pgs_filt_cis, idcol = "eid")

pgs_filt_trans= fread("data/FG_filtered_trans_pPGS/trans_genoresiduals.tsv") 
###scale 
pgs_filt_trans = scale_df(pgs_filt_trans,  idcol = "eid")

residuals_filt_gw= fread("data/FG_filtered_gw_pPGS/gw_residuals_output_nocov.tsv") %>% select(!(AggAbsResid)) 
###scale 
residuals_filt_gw = scale_df(residuals_filt_gw,  idcol = "IndID")

pgs_filt_gw= fread("data/FG_filtered_gw_pPGS/gw_genoresiduals.tsv") 
###scale 
pgs_filt_gw = scale_df(pgs_filt_gw,  idcol = "eid")

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

### INCIDENT LOGISTIC REGRESSION #################
incident_logif_pairs <- function(prot_df, idcol,
                                 endpoint_date_suffix = "_DATE",
                                 sep = "__") {
  
  stopifnot(idcol %in% names(prot_df))
  
  # available disease date columns in intervene_fi
  date_cols <- names(intervene_fi %>% select(contains(endpoint_date_suffix)))
  
  res_all <- list()
  
  ## loop over colums with protein-disease pairs
  for (score_col in setdiff(names(prot_df), idcol)) {
    
    # Expect "PROTEIN__ENDPOINT"
    parts <- str_split_fixed(score_col, fixed(sep), 2)
    prot  <- parts[1]
    endp  <- parts[2]
    
    # map endpoint -> date column name 
    dis_date <- paste0(endp, endpoint_date_suffix)
    
    df_base <- intervene_fi %>%
      select(ID, END_OF_FOLLOWUP, baseline, age, sex, all_of(dis_date))
    
    p <- prot_df %>%
      select(ID = all_of(idcol), score = all_of(score_col))
    
    df_merge <- merge(df_base, p, by = "ID")
    
    # follow-up time 
    df_merge <- df_merge %>%
      mutate(fol = difftime(.data[[dis_date]], baseline, units = "weeks")/52.25) ### note that I converted in years!!
    
    # excluding all prevalent cases or incident cases recorded within the first 6 months of follow-up
    df_filt <- df_merge %>% filter(fol > 0.5 | is.na(fol))
    
    ### binary outcome, 1= incidence, 0 = right censored
    df_bin <- df_filt %>%
      mutate(y = ifelse(is.na(.data[[dis_date]]), 0, 1))
    
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

### results CIS-RESIDUALS FILTERED 
results_cis_filt = incident_logif_pairs(prot_df = residuals_filt_cis, idcol = "IndID")
results_cis_filt$protein= gsub("StandResid_","", results_cis_filt$protein)
results_cis_filt$Model= "cis-residualsFILT"


### results TRANS-RESIDUALS FILTERED 
results_trans_filt = incident_logif_pairs(prot_df = residuals_filt_trans, idcol = "IndID")
results_trans_filt$protein= gsub("StandResid_","", results_trans_filt$protein)
results_trans_filt$Model= "trans-residualsFILT"

### results CIS-PGS FILTERED 
results_cisPGS_filt = incident_logif_pairs(prot_df = pgs_filt_cis, idcol = "eid")
results_cisPGS_filt $Model= "cis-PGS_FILT"

### results TRANS-PGS FILTERED 
results_transPGS_filt = incident_logif_pairs(prot_df = pgs_filt_trans, idcol = "eid")
results_transPGS_filt$Model= "trans-PGS_FILT"

### results GENOME WIDE RESIDUALS FILTERED 
results_filt = incident_logif_pairs(prot_df = residuals_filt_gw, idcol = "IndID")
results_filt$Model= "residuals_FILT"

### results genomewide-PGS FILTERED 
results_PGS_filt = incident_logif_pairs(prot_df = pgs_filt_gw, idcol = "eid")
results_PGS_filt$Model= "PGS_FILT"

results_all= fread("data/cistransresiduals_logisticregressions.tsv") 

### merge all 
all= bind_rows(results_filt, results_cis_filt, results_trans_filt,
               results_PGS_filt, results_cisPGS_filt, results_transPGS_filt,
               results_all)

all = all %>% rename(zscore = `z value`,
                     pval= "Pr(>|z|)",
                     se= "Std. Error")
all$protein= gsub("StandResid_","",all$protein)


write.table(all, "data/V2_revision_logisticmodels_protdisease.tsv", quote=F, sep="\t", row.names = F)







