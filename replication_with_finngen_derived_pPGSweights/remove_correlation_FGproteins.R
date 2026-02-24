library(data.table)
library(dplyr)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

ourprot= fread("data/R2_94proteins_olinksoma.tsv")
ukbolink= fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/omics_phenotypes/olink_data.txt")
prot_tr= fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/omics_phenotypes/ukb_olink_proteomics_translations.tsv")
our_id= fread("data/idefix_olinksoma/pheno_residuals.tsv")

ukb_olink= ukbolink %>% filter(ins_index==0)
ukb_olink= merge(ukb_olink, prot_tr, by="protein_id")

ukb_94= ukb_olink %>%
  filter(gene_name %in% ourprot$protname & eid %in% our_id$eid)### just check if there is >1 panel
ukb_94_mat= dcast(ukb_94, eid ~ gene_name, value.var = "result")

fg_prot= fread("/scratch/project_2007428/users/Zhiyu/TmpProtein/KeepProt_h2Prune",header=F)
ukb_172= ukb_olink %>%
  filter(gene_name %in% fg_prot$V1 & eid %in% our_id$eid)

ukb_fgprot= ukb_172 %>%
  filter(!(gene_name %in%
             intersect(colnames(ukb_94_mat),unique(ukb_172$gene_name))))

ukb_fgprot_mat= dcast(ukb_fgprot, eid ~ gene_name, value.var = "result")

c= cor(ukb_fgprot_mat[,-"eid"], ukb_94_mat[,-"eid"],
       use = "pairwise.complete.obs", method = "pearson")

high_corr <- which(abs(c) > 0.5, arr.ind = TRUE)
rownames(c)[high_corr[, "row"]]   # rownames (from FG proteins) --> CXCL5 is present twice

final_list= rownames(c)[!(rownames(c) %in% unique(rownames(high_corr)))] ## remove correlated proteins with our 94list
write.table(data.frame(final_list), file="data/nocor_finallist_FGprot.txt",
            quote=F,row.names = F,sep="\t")
