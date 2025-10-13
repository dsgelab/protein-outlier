library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggnewscale)  # allows multiple fill scales
library(plotly)
library(stringr)
library(nnet)
library(foreach)
library(doParallel)
library(fastDummies)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")


envexp= fread("/scratch/project_2007428/users/Zhiyu/Pheno/PhenoOut/EnvPheno")

resid= fread("data/idefix_olinksoma/residuals_output_nocov.tsv") ## LOAD RESIDUALS NO COV
setDT(resid[, (names(resid)[-1]) := lapply(.SD, scale), .SDcols = names(resid)[-1]])  ## scale residuals 
resid$AggAbsResid=NULL

base_pheno= fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv")
base_pheno= base_pheno %>% select(eid,testing_site= `54-0.0`, BMI= `21001-0.0`, age=`21022-0.0`, sex= `31-0.0`)

bloodcov=fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_blood_biochem_phenotypes.tsv")
bloodcov= bloodcov %>% select(eid, platelet_count= '30080-0.0', creatinine="30700-0.0") 

envexp = envexp %>% filter(eid %in% resid$IndID)

nacount= colSums(is.na(envexp))

envs_common= names(nacount[nacount<1000]) 

description= fread("/scratch/project_2007428/users/Zhiyu/Pheno/PhenoList/PhenoMemo",fill=TRUE,header = FALSE)
description= description %>% filter(!V2=="Household") %>%
  mutate(Category = ifelse(grepl("^#", V1), V2, NA)) %>%
  fill(Category, .direction = "down") %>%
  rename(code=V1,descr=V2)
description= description[!grepl("^#", description[[1]])]

description= description %>% filter(code %in% gsub("-0.0","",envs_common))

env= envexp %>% select(eid, envs_common)
colnames(env)= gsub("-0.0","",colnames(env))
#colnames(env)=  description$descr[match(names(env), description$code, nomatch="eid")]

#### see the distibution of phenotypes to separate categorical from ordinal ones. 
dist=data.frame(
  exp = colnames(env)[-1],
  N = sapply(env[,-1], uniqueN)
)


df= merge(env, resid, by.x="eid", by.y="IndID")
### add other variables 
df= merge(base_pheno, df, by="eid")
df= merge(bloodcov, df, by="eid")


####### RUN REGRESSION, CONVERT IN DUMMY VARIABLES
results= data.frame()
cl = makeCluster(20)
registerDoParallel(cl)
res= foreach(prot = colnames(df)[93:ncol(df)],
             .combine=bind_rows, .packages=c("nnet","broom","dplyr","data.table","fastDummies")) %:%
  foreach(e = colnames(df)[c(2:5,8:92)], .combine=bind_rows) %dopar% {
    
    df_filt= df
    
    if (uniqueN(df[[e]])<20 | e=="testing_site") {
      
      tab = table(df_filt[[e]])
      valid_levels = names(tab[tab >= 100])
      
      df_filt = df_filt %>% filter(.data[[e]] %in% valid_levels) ### filter out levels with < 100 entries
      
      #drop unused factor levels if needed
      if (is.factor(df_filt[[e]])) {
        df_filt[[e]] = droplevels(df_filt[[e]])
      }
      
      # Use dummy_cols to one-hot encode the column `e`
      df_filt= dummy_cols(
        df_filt,
        select_columns = e,
        remove_selected_columns = T,
        remove_first_dummy = F,
        ignore_na = T ### does not create a dummy for NAs
      )
      
      # Get the names of the created dummy columns
      dummy_cols_created = grep(paste0("^", e, "_"), names(df_filt), value = T)
      
      # Run logistic regression for each dummy column as response
      results = lapply(dummy_cols_created, function(dummy_env) {
        formula= as.formula(paste0("`",dummy_env,"` ~", prot, "+ age + sex"))
        model= glm(formula, data = df_filt, family = "binomial")
        tidy(model) %>% filter(!term %in% c("(Intercept)","age","sex")) %>%
          mutate(class = dummy_env, exposure = e)
      })
      
      bind_rows(results)
      
    } else {
      # For continuous exposures, run linear regression
      form = paste0("scale(`",e,"`) ~ ",prot," + age + sex")
      m= lm(as.formula(form), data=df)
      results= tidy(m) %>% filter(!term %in% c("(Intercept)","age","sex")) %>%
        mutate(class = NA, exposure = e)
    }
    
    bind_rows(results)
  }
stopCluster(cl)
#write.table(res, file="data/environment_residuals_dummy_reg.tsv",quote=F,row.names=F,sep="\t")


############################## with proteins 
prot= fread("data/idefix_olinksoma/pheno_residuals.tsv")
#scale
setDT(prot[, (names(prot)[-1]) := lapply(.SD, scale), .SDcols = names(prot)[-1]])

df2= merge(env, prot, by="eid")
### add other variables 
df2= merge(base_pheno, df2, by="eid")
df2= merge(bloodcov, df2, by="eid")

res_prot= foreach(prot = colnames(df2)[93:ncol(df2)],
             .combine=bind_rows, .packages=c("nnet","broom","dplyr","data.table","fastDummies")) %:%
  foreach(e = colnames(df2)[c(2:5,8:92)], .combine=bind_rows) %dopar% {
    
    df_filt= df2
    
    if (uniqueN(df2[[e]])<20 | e=="testing_site") {
      
      tab = table(df_filt[[e]])
      valid_levels = names(tab[tab >= 100])
      
      df_filt = df_filt %>% filter(.data[[e]] %in% valid_levels) ### filter out levels with < 100 entries
      
      #drop unused factor levels if needed
      if (is.factor(df_filt[[e]])) {
        df_filt[[e]] = droplevels(df_filt[[e]])
      }
      
      # Use dummy_cols to one-hot encode the column `e`
      df_filt= dummy_cols(
        df_filt,
        select_columns = e,
        remove_selected_columns = T,
        remove_first_dummy = F,
        ignore_na = T ### does not create a dummy for NAs
      )
      
      # Get the names of the created dummy columns
      dummy_cols_created = grep(paste0("^", e, "_"), names(df_filt), value = T)
      
      # Run logistic regression for each dummy column as response
      results = lapply(dummy_cols_created, function(dummy_env) {
        formula= as.formula(paste0("`",dummy_env,"` ~", prot, "+ age + sex"))
        model= glm(formula, data = df_filt, family = "binomial")
        tidy(model) %>% filter(!term %in% c("(Intercept)","age","sex")) %>%
          mutate(class = dummy_env, exposure = e)
      })
      
      bind_rows(results)
      
    } else {
      # For continuous exposures, run linear regression
      form = paste0("scale(`",e,"`) ~ ",prot," + age + sex")
      m= lm(as.formula(form), data=df2)
      results= tidy(m) %>% filter(!term %in% c("(Intercept)","age","sex")) %>%
        mutate(class = NA, exposure = e)
    }
    
    bind_rows(results)
  }
stopCluster(cl)
#write.table(res_prot, file="data/environment_proteins_dummy_reg.tsv",quote=F,row.names=F,sep="\t")


#### rerun regression by selecting most significant results for one example. 
#protein= CCL15, envs= current tobacco smoking (1239), alcohol intake frequency (1558) ,
#                      Attendance/disability/mobility allowance(6146) ,
#                      Mayor dietary changes in the last 5 years (1538),
#                      Own or rent accomodation lived in (680) 
### residuals
df_f= df %>% select(eid,age,sex,
                    smoking=`1239`,
                    alcohol=`1558`,
                    allowance=`6146`,
                    dietary_change=`1538`,
                    accomodation= `680`,
                    StandResid_CCL15)
prot= fread("data/idefix_olinksoma/pheno_residuals.tsv") %>% select(eid,CCL15)
#scale
setDT(prot[, (names(prot)[-1]) := lapply(.SD, scale), .SDcols = names(prot)[-1]])

df_f= merge(df_f,prot,by="eid")

######## CKD 
### file with first occurrences of INTERVENE phenotypes 
intervene_fi= arrow::read_parquet("/scratch/project_2007428/projects/prj_010_phenotype_file/share/ukb78537_INTERVENE_phenotype-2024-08-01.parquet")
intervene_fi= intervene_fi %>% select(ID,SEX,END_OF_FOLLOWUP, N14_CHRONKIDNEYDIS_DATE)
### filter for our individuals
intervene_fi = intervene_fi %>% filter(ID %in% df_f$eid)
### filter phenotypes for the intervene. note that from the 39 endpoints, BMI and Covid-19 severity are not present as first occurrences, so the final number is 37 
definitions= fread("data/Intervene_flagship_endpoint_collection_Definitions.tsv") 
definitions = definitions %>% select(Endpoint,`FinnGen endpoint`)
base_pheno= fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv") %>%
  select(ID=eid, baseline='53-0.0', death='40000-0.0')
intervene_fi= merge(intervene_fi,base_pheno,by="ID")

df_merge = merge(df_f, intervene_fi,by.y="ID",by.x="eid")

df_merge= df_merge %>% mutate(fol = difftime(N14_CHRONKIDNEYDIS_DATE, baseline, units = "weeks")/52.25) ###note that I converted in years!!
## excluding all prevalent cases or incident cases recorded within the first 6 months of follow-up
df_filt= df_merge %>% filter(fol>0.5 | is.na(fol))
### binary outcome, 1= incidence, 0 = right censored
ckd_incident = df_filt %>% mutate(CKD := ifelse(is.na(N14_CHRONKIDNEYDIS_DATE),0,1))


# Recode enviromental exposure functions 
recode_smoking <- function(x) {           ######## Current tobacco smoking
  case_when(
    x %in% c(1, 2) ~ 1,
    x == 0 ~ 0,
    x == -3 ~ NA_real_
  )
}
recode_alcohol <- function(x) {     ###### Regular drinker
  case_when(
    x %in% c(1,2,3) ~ 1,
    x %in% c(4,5,6) ~ 0,
    x == -3 ~ NA_real_
  )
}
recode_allowance <- function(x) {                 ############## disability allowance
  case_when(
    x %in% c(1, 2, 3) ~ 1,
    x == -7 ~ 0,
    x %in% c(-1, -3) ~ NA_real_
  )
}
recode_dietary <- function(x) {          ##### major dietary changes BECAUSE OF ILLNESS
  case_when(
    x %in% c(0, 2) ~ 0,
    x == 1 ~ 1,
    x == -3 ~ NA_real_
  )
}

recode_accomodation <- function(x) {              ##########Public social house rental
  case_when(
    x %in% c(1, 2, 4, 5, 6, -7) ~ 0,
    x == 3 ~ 1,
    x == -3 ~ NA_real_
  )
}

# List of exposures with their recode functions and names
exposures= list(
  list(var = "smoking", name = "Current tobacco smoking", recode = recode_smoking),
  list(var = "alcohol", name = "Regular alcohol drinker", recode = recode_alcohol),
  list(var = "allowance", name = "Disability allowance", recode = recode_allowance),
  list(var = "dietary_change", name = "Major dietary changes in the last 5y because of illness", recode = recode_dietary),
  list(var = "accomodation", name = "Public/social house rental", recode = recode_accomodation)
)

# Predictors to run regression with
predictors= c("StandResid_CCL15", "CCL15")

library(purrr)
# Run all analyses and combine
results = map_dfr(exposures, function(e) {
  map_dfr(predictors, function(pred) {
    run_analysis(df_f, e$var, e$name, e$recode, pred)
  })
})

#### add the CKD - Tobacco smoking regression
results= rbind(results,
               run_analysis(ckd_incident, "CKD", "CKD", recode_smoking, "smoking")
)

results$model = ifelse(results$model=="StandResid_CCL15","genetically adjusted protein",
                       ifelse(results$model=="CCL15","unadjusted protein","none"))
write.table(results, file="data/sigenv_reg_CCL15.tsv",quote=F,row.names=F,sep="\t")



