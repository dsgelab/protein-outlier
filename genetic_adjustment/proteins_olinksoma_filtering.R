library(data.table)
library(reshape2)
library(dplyr)
library(caret)
args=commandArgs(T)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

r2olink= args[1]
r2soma= args[2]
r2_cutoff= args[3]
r_cutoff= args[4]
prot_output= args[5]
id_output= args[6]
pheno_output= args[7]

#r2_olink=fread("data/R2_pPRSprot_nocov.tsv")
r2_olink= fread(r2olink)
olink_selprot= fread("data/omicspred_selected_proteins.tsv")
r2_olink= merge(r2_olink,olink_selprot, by.x="protname", by.y="Gene")
r2_olink= r2_olink %>% dplyr::rename('omicspred_id'='OMICSPRED ID') 
r2_olink$type= "olink"

#r2_soma= fread("data/Somalogic_R2_pPRSprot_nocov.tsv")
r2_soma= fread(r2soma)
soma_selprot= fread("data/somalogic_selected_proteins.csv") %>% dplyr::select(Gene, 'UniProt ID') 
r2_soma= merge(r2_soma, soma_selprot, by.x="protname", by.y="Gene")
r2_soma= r2_soma[!duplicated(r2_soma),]
r2_soma$type= "somalogic"

all_r2= rbind(r2_olink,r2_soma)
filt_r2= all_r2[all_r2$R2>as.numeric(r2_cutoff),] ### filter for R2 > cutoff
### here I am taking the protein with best prediction (higher R2), since I have duplicated ones from olink and somalogic
filt_r2= filt_r2 %>% group_by(protname) %>% slice(which.max(R2)) %>% ungroup()

correl= function(cutoff) {
  dir= "/scratch/project_2007428/projects/prj_100_pprs_discordance/data/"
  files=list.files(c(paste0(dir,"pred_real_protvalues"),paste0(dir,"Somalogic_pred_real_protvalues")), 
                   full.names = T)
  keywords= paste0(filt_r2$protname,"_",filt_r2$omicspred_id,"_predrealprot.tsv")
  sel_files= files[sapply(files, function(x) any(sapply(keywords, function(y) grepl(y, x))))]
  
  alldata= data.frame()
  for (f in sel_files) {
    data= read.table(f, header=F)
    colnames(data)= c("sample","predicted","real")
    data$protein= strsplit(strsplit(f,"/")[[1]][8], "_OPG")[[1]][1]
    alldata=rbind(data[c("sample","real", "protein")],alldata)
  }
  df=reshape2::dcast(alldata, sample ~ protein, value.var = "real", drop=T)
  df=df[,-1]
  cor_matrix=cor(df, use = "pairwise.complete.obs")
  
  a=caret::findCorrelation(cor_matrix, cutoff=cutoff, names=T, verbose = TRUE )
  nocor= filt_r2[!(filt_r2$protname %in% a),]
  return(nocor)
}

allprot= correl(as.numeric(r_cutoff))

prot_omicspred= allprot %>% select(protname, omicspred_id)
write.table(prot_omicspred, file= prot_output, quote=F, sep="\t", row.names = F)

imp=fread("/scratch/project_2007428/projects/prj_004_sex_prediction/data/processed/proteomics_imputed_2024-02-06.csv")
prot_tr= fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/omics_phenotypes/ukb_olink_proteomics_translations.tsv")
agesex= fread("/scratch/project_2007428/projects/prj_100_pprs_discordance/data/EURprotindv_agesex.tsv")
colnames(agesex)=c("sample","sex","age")
h= merge(allprot,prot_tr, by.x="protname", by.y="gene_name")
imp_sel= imp %>% select(c(as.character(h$protein_id),"eid")) ###select proteins
imp_sel= imp_sel %>% filter(eid %in% agesex$sample) ### select individuals
row.names(imp_sel)= imp_sel$eid
pheno= reshape2::melt(imp_sel, id="eid")
pheno= merge(pheno, prot_tr, by.x="variable", by.y="protein_id")
pheno= pheno %>% select(eid,value,gene_name)

###remove unmatched genetic sex - sex 
matched_sex=fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv") 
matched_sex= matched_sex %>% select(eid, `22001-0.0`, `31-0.0`)
colnames(matched_sex)= c("eid", "genetic_sex", "sex")
matched_sex= matched_sex %>% filter(eid %in% unique(pheno$eid))
matched_sex= matched_sex[matched_sex$genetic_sex==matched_sex$sex,] ### take only if genetic sex==sex
pheno= pheno %>% filter(eid %in% matched_sex$eid)

### add age and and sex
pheno= merge(pheno, agesex, by.x="eid",by.y="sample")
indv= pheno %>% select(eid,sex,age) 
indv= indv[!(duplicated(indv)),]
write.table(indv, file= id_output, quote=F, sep="\t", row.names = F)

pheno[pheno$sex==0,]$sex= "Female"
pheno[pheno$sex==1,]$sex= "Male"
colnames(pheno)= c("ID","VALUE","TRAIT","SEX","AGE")

write.table(pheno, file= pheno_output, quote=F, sep="\t", row.names = F)


