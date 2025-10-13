library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")
args=commandArgs(T)
indv= args[1]
proteins= args[2]
out_proteins= args[3]
out_pprs= args[4]

#### PHENO file
### individuals
indv= fread(indv) %>% select(ID)
### proteins
prot= fread(proteins)

### protein code and protein name 
prot_tr= fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/omics_phenotypes/ukb_olink_proteomics_translations.tsv")
prot_tr= prot_tr %>% filter(gene_name %in% prot$protname)

imp_prot=fread("/scratch/project_2007428/projects/prj_004_sex_prediction/data/processed/proteomics_imputed_2024-02-06.csv")
imp_prot= imp_prot %>% filter(eid %in% unique(indv$ID)) 

selimp_prot <- imp_prot %>%
  select(c("eid",as.character(prot_tr %>% pull(protein_id))))

name_mapping <- setNames(prot_tr$gene_name, prot_tr$protein_id)
selimp_prot <- selimp_prot %>%
  rename_with(~ ifelse(. %in% names(name_mapping), name_mapping[.], .), .cols = 2:ncol(selimp_prot))
write.table(selimp_prot, file=out_proteins, sep="\t",
            quote=F,col.names = T,row.names = F)

### pgs file
pprs= list.files(c("data/pPRS","data/pPRS_Somalogic"), pattern = ".sscore", full.names = T)
prot_= paste0(prot$protname,"_",prot$omicspred_id)
filt_pprs=pprs[sapply(pprs, function(pprs) any(sapply(prot_, function(pattern) grepl(pattern, pprs))))]
allpprs= data.frame(eid= unique(indv$ID))
for (file in filt_pprs) {
  df= fread(file) %>% select(IID,SCORE1_AVG)
  prot= str_extract(file, "(?<=EUR_)[^_]+(?=_OPG)")
  colnames(df)= c("eid",prot)
  df= df %>% filter(eid %in% unique(indv$ID))
  allpprs= merge(allpprs,df, by="eid")
}
write.table(allpprs, file=out_pprs, sep="\t",
            quote=F,col.names = T,row.names = F)



