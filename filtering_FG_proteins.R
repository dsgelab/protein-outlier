library(data.table)
library(dplyr)
library(tidyr)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

files=list.files("data/FG_pred_real_protvalues/")

alldata= data.frame()
r2tab= data.frame()
for (f in files) {
  r2table=data.frame(ncol=7)
  r2table$protname= strsplit(f, "_")[[1]][1]
  data= fread(paste0("data/FG_pred_real_protvalues/",f), header = F)
  colnames(data)= c("sample","predicted","real")
  model=lm(real ~ predicted, data = data)
  r2table$sample= nrow(data)
  r2table$R2=summary(model)$r.squared
  r2table$adj_R2=summary(model)$adj.r.squared
  r2table$beta= summary(model)$coefficients[2, "Estimate"]
  r2table$se= summary(model)$coefficients[2, "Std. Error"]
  r2table$pvalue= summary(model)$coefficients[2, "Pr(>|t|)"]
  r2tab=rbind(r2tab, r2table)
  
}
r2tab$ncol=NULL

write.table(r2tab, file="data/FGproteins_R2.tsv",quote=F,sep="\t",row.names = F)

####### load
r2= fread("data/FGproteins_R2.tsv")
r2_filt= r2 %>% filter(R2>0.2)

########## PHENO FOR RESIDUALS  
### protein code and protein name 
prot_tr= fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/omics_phenotypes/ukb_olink_proteomics_translations.tsv")
prot_tr= prot_tr %>% filter(gene_name %in% r2_filt$protname)

imp_prot=fread("/scratch/project_2007428/projects/prj_004_sex_prediction/data/processed/proteomics_imputed_2024-02-06.csv")
#imp_prot= imp_prot %>% filter(eid %in% unique(indv$ID)) 

selimp_prot <- imp_prot %>%
  select(c("eid",as.character(prot_tr %>% pull(protein_id))))

name_mapping <- setNames(prot_tr$gene_name, prot_tr$protein_id)
selimp_prot <- selimp_prot %>%
  rename_with(~ ifelse(. %in% names(name_mapping), name_mapping[.], .), .cols = 2:ncol(selimp_prot))

ids= fread("data/indv_agesex_olinksoma.tsv")

selimp_prot= selimp_prot %>% filter(eid %in% ids$eid)

write.table(selimp_prot, file="data/FG_pheno_residuals.tsv", sep="\t",
            quote=F,col.names = T,row.names = F)


######### GENO FOR RESIDUALS
### pgs file
pprs= list.files("data/pPRS_FGproteins", pattern = ".sscore", full.names = T)
prot_= r2_filt$protname
filt_pprs=pprs[sapply(pprs, function(pprs) any(sapply(prot_, function(pattern) grepl(pattern, pprs))))]
allpprs= data.frame(eid= selimp_prot$eid)
for (file in filt_pprs) {
  df= fread(file) %>% select(IID,SCORE1_AVG)
  prot= stringr:: str_extract(file, "(?<=EUR_).*?(?=\\.sscore)")
  colnames(df)= c("eid",prot)
  df= df %>% filter(eid %in% selimp_prot$eid)
  allpprs= merge(allpprs,df, by="eid")
}

write.table(allpprs, file="data/FG_geno_residuals.tsv", sep="\t",
            quote=F,col.names = T,row.names = F)





