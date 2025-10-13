library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(patchwork)
library(ggrepel)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

allproteins=fread("data/FG_pheno_residuals.tsv")
###scale proteins
setDT(allproteins[, (names(allproteins)[-1]) := lapply(.SD, scale), .SDcols = names(allproteins)[-1]])

allresid= fread("data/FG_residuals_output_nocov.tsv")
###scale residuals
setDT(allresid[, (names(allresid)[-1]) := lapply(.SD, scale), .SDcols = names(allresid)[-1]])

allprs= fread("data/FG_geno_residuals.tsv")
### scale PRS
setDT(allprs[, (names(allprs)[-1]) := lapply(.SD, scale), .SDcols = names(allprs)[-1]])

### file with first occurrences of INTERVENE phenotypes 
intervene_fi= arrow::read_parquet("/scratch/project_2007428/projects/prj_010_phenotype_file/share/ukb78537_INTERVENE_phenotype-2024-08-01.parquet")
intervene_fi= intervene_fi %>% select(ID,SEX,END_OF_FOLLOWUP, contains("_DATE"))
### filter for our individuals
intervene_fi = intervene_fi %>% filter(ID %in% allproteins$eid)

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

################################# INCIDENT LOGISTIC REGRESSION
incident_logif= function(prot_df, idcol) {
  allendpoints = data.frame()
  for (dis in colnames(intervene_fi %>% select(contains("_DATE")))) {
    df= intervene_fi %>% select(ID, END_OF_FOLLOWUP, dis, baseline, age, sex) 
    allresults= data.frame()
    for (prot in colnames(prot_df)[-1]) {
      p= prot_df %>% select(ID= idcol, prot)
      df_merge= merge(df, p, by="ID")
      df_merge= df_merge %>% mutate(fol = difftime(!!as.name(dis), baseline, units = "weeks")/52.25) ###note that I converted in years!!
      
      ## excluding all prevalent cases or incident cases recorded within the first 6 months of follow-up
      df_filt= df_merge %>% filter(fol>0.5 | is.na(fol))
      
      ### binary outcome, 1= incidence, 0 = right censored
      df_bin= df_filt %>% mutate(!!as.name(dis) := ifelse(is.na(!!as.name(dis)),0,1))
      
      ### Logistic regression
      formula= as.formula(paste0(dis,"~ `",prot,"` + age + sex"))
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


incident_prot= incident_logif(allproteins,"eid") 
incident_prot= incident_prot %>% rename(estimate_prot= Estimate, SE_prot= `Std. Error`,
                                        zvalue_prot= `z value`, pval_prot =`Pr(>|z|)`)
write.table(incident_prot, file="data/FG_incident_logreg_singleproteins.tsv", sep="\t", quote=F,
            row.names = F, col.names = T)

incident_resid= incident_logif(allresid, "IndID")
incident_resid = incident_resid %>% rename(estimate_resid= Estimate, SE_resid= `Std. Error`,
                                           zvalue_resid= `z value`, pval_resid =`Pr(>|z|)`)
write.table(incident_resid, file="data/FG_incident_logreg_singleresiduals.tsv", sep="\t", quote=F,
            row.names = F, col.names = T)

incident_prs= incident_logif(allprs, "eid")
incident_prs= incident_prs %>% rename(estimate_PRS= Estimate, SE_PRS= `Std. Error`,
                                      zvalue_PRS= `z value`, pval_PRS =`Pr(>|z|)`)
write.table(incident_prs, file="data/FG_incident_logreg_singlePRS.tsv", sep="\t", quote=F,
            row.names = F, col.names = T)
#-----------------------------------------------------------------------------------------------------------


