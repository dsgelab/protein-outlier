library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(ggrepel)
library(ggnewscale)  # allows multiple fill scales
library(htmlwidgets)
library(plotly)
library(stringr)
library(nnet)
library(foreach)
library(doParallel)
library(fastDummies)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")


### here I'm running exposome regressions just using the cis-adjusted proteins

envexp= fread("/scratch/project_2007428/users/Zhiyu/Pheno/PhenoOut/EnvPheno")

### cis residuals
resid= fread("data/cis_residuals_output_nocov.tsv") ## LOAD RESIDUALS NO COV
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
write.table(res, file="data/environment_cisresiduals_dummy_reg.tsv",quote=F,row.names=F,sep="\t")


