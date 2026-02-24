library(dplyr)
library(data.table)
library(argparse)

# Set a safe temporary directory for bigstatsr file-backed matrices
tmpdir <- Sys.getenv("TMPDIR", unset = "/tmp")
if (!dir.exists(tmpdir)) dir.create(tmpdir, recursive = TRUE)
options(bigstatsr.backingpath = tmpdir)
cat("Backing files will be written to:", getOption("bigstatsr.backingpath"), "\n")

setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")
source("scripts/LDpredest_per_chr_func.R")
parser= ArgumentParser(description = "LDpred weights for PGS")

parser$add_argument("--prot", required = T)
parser$add_argument("--outfile", required = T)
parser$add_argument("--h2", required = T)

args = parser$parse_args()

## SumHEM per window
df_map <- readRDS("data/LDref_LDpred/map_hm3_plus.rds")
df_map$MAF <- ifelse(df_map$af_UKBB > 0.5,1 - df_map$af_UKBB,df_map$af_UKBB)
df_map_unduplicated <- distinct(df_map,rsid,.keep_all = TRUE)
  
print(args$prot)
## reading gwas
  
df_gwas <- fread(paste0("/scratch/project_2007428/data/summary_stats/FinnGenProtSelect/",args$prot,".txt.gz"))
head(df_gwas)

### convert SNPs in rsid
snps= fread("data/FinnGenSNP")
snps= snps %>% mutate(SNP= paste0("chr",gsub(":","_",SNP)))

##add the rsid column from snps file 
df_gwas= merge(df_gwas, snps, by.x="ID", by.y="SNP")

df_gwas1 <- df_gwas %>% 
      rename("chr"="CHR", "pos_hg38"="POS", "beta"="BETA", "se"="SE", 
             "OA" = "REF", "EA" = "ALT", "EAF" = "ALT_FREQ") %>%
  select(rsid,chr,pos_hg38,OA,EA,EAF,beta,se,N) %>%
  distinct(rsid,.keep_all = TRUE)
  
  ## match allele
  
df_gwas1 <- df_gwas1 %>% inner_join(df_map_unduplicated %>% select(rsid,chr,pos,pos_hg38,a0,a1,af_UKBB))

with(df_gwas1,table(OA == a0))

with(df_gwas1,table(OA == a0 & EA == a1))
  
df_gwas1 <- df_gwas1[EA == a1,]
  
with(df_gwas1,table(OA == a0 & EA == a1))
  
head(df_gwas1)
## LDpred
  
LDpred_res <- LDpredest(df_gwas = df_gwas1, 
                        
                        df_map = df_map %>% select(rsid,chr,pos,ld,MAF), 
                        
                        MAF_threshold = 0.01,
                      
                        LD_path="data/LDref_LDpred/ldref_hm3_plus/", 
                        
                        h2input= as.numeric(args$h2),
                        
                        NCORES = 8)
  
  
  
LDpred_res <- merge(LDpred_res,df_gwas1[,c("rsid","chr","pos","a0","a1")])
  
setcolorder(LDpred_res, c(
    
  "rsid", "chr", "pos", "a0", "a1", "beta",
    
  setdiff(colnames(LDpred_res), c("rsid", "chr", "pos", "a0", "a1", "beta"))
    
))

  
setorder(LDpred_res,chr,pos)
  
  # original LDpred results
saveRDS(LDpred_res, file = paste0("data/LDpred_res/df_est.LDpred.MAF01.",args$prot,".rds"))
  
LDpred_res[,b_auto_LDpred_not_scale := b_auto_LDpred * scale]
LDpred_res_subset <- LDpred_res[,.(rsid,a1,b_auto_LDpred_not_scale)]
  
# output SNP weights
fwrite(LDpred_res_subset,file = args$outfile,
       sep = "\t",col.names = F)


