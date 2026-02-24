library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(ggrepel)
library(biomaRt)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

r2all= fread("data/R2_94proteins_olinksoma.tsv")
r2all = r2all %>% dplyr::select(protein= protname,R2,type)

#### description of the endpoints
definitions= fread("data/Intervene_flagship_endpoint_collection_Definitions.tsv") 
definitions = definitions %>% dplyr::select(Endpoint,`FinnGen endpoint`)

#### LOAD INCIDENT LOGISTIC REGRESSIONS
inc_prot= fread("data/incident_logreg_singleproteins.tsv")
inc_resid= fread("data/incident_logreg_singleresidualsnocov.tsv") %>% filter(protein != "AggAbsResid") %>% mutate(protein= gsub("StandResid_","",protein))
inc_prs=  fread("data/incident_logreg_singlePRS.tsv")

all_incident= inc_prot %>% inner_join(inc_resid, by=c("endpoint","protein")) %>% inner_join(inc_prs, by=c("endpoint","protein")) 
all_incident= all_incident %>% mutate(endpoint= gsub("_DATE","",endpoint))

### filter phenotypes for the intervene. note that from the 39 endpoints, BMI and Covid-19 severity are not present as first occurrences, so the final number is 37
all_incident= merge(all_incident,definitions, by.x="endpoint",by.y="FinnGen endpoint")
all_incident= merge(all_incident,r2all,by="protein")
#all_incident_sig= all_incident %>% filter(pval_prot<0.05 | pval_resid<0.05)
all_incident= all_incident %>% mutate(p_adjprot= p.adjust(pval_prot, method="bonferroni"),
                                      p_adjresid= p.adjust(pval_resid, method="bonferroni"))
#### beta difference ( residuals - protein )
all_incident = all_incident %>% mutate(diffresprot= estimate_resid-estimate_prot)
### SE of the difference 
all_incident = all_incident %>% mutate(se_diffresprot= sqrt(SE_resid^2 + SE_prot^2))
### z-score of the difference
all_incident  = all_incident  %>% mutate(z_diffresprot = diffresprot/se_diffresprot)
### two-tailed p-value of the z-score 
all_incident  = all_incident  %>% mutate(p_diffresprot= 2 * (1 - pnorm(abs(z_diffresprot))))

all_incident= all_incident %>% mutate(sig_bonf= ifelse(p_adjprot < 0.05 | p_adjresid < 0.05,"yes","no"))

#### significant pairs after bonferroni correction
all_incident_padj= all_incident %>% filter(p_adjprot < 0.05 | p_adjresid < 0.05)

##### pairs significant at z-test, keeping outliers 
moresig= all_incident_padj[all_incident_padj$p_diffresprot<0.05 & abs(all_incident_padj$z_diffresprot)>3,] 
#proteins= unique(moresig$protein)
proteins= unique(all_incident$protein)

# Connect to the GRCh37 archive Ensembl server
ensembl37=useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl",host="https://grch37.ensembl.org")

# Query chromosome and position on GRCh37
proteins_pos=getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = proteins,
  mart = ensembl37
)
proteins_pos = proteins_pos %>% filter(!(grepl("H", chromosome_name)))

olink= fread("data/Olink_trait_validation_results_with_OMICSPRED_ID.csv")
olink= olink %>% dplyr::select(omicspred_id=`OMICSPRED ID`,protein= Gene) %>% filter(protein %in% proteins)
olink$type= "olink"

somalogic= fread("data/Somalogic_R2_pPRSprot_nocov.tsv") %>% dplyr::select(protein=protname,omicspred_id) %>% filter(protein %in% proteins)
somalogic$type= "somalogic"

olink_soma= rbind(olink,somalogic)
olink_soma=merge(olink_soma,r2all,by=c("protein","type"))

protein_pos= merge(proteins_pos,olink_soma,by.x="hgnc_symbol",by.y="protein") %>% rename(protein=hgnc_symbol)

write.table(protein_pos, file="data/iLR_outliers_prot_pos.tsv",quote=F,sep="\t",row.names=F)
