library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

pairs <- fread("data/protein_loci_3478pairs.tsv")

count_snps <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  length(count.fields(path)) - 1L  # subtract header
}

count_snps_nohead <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  length(count.fields(path)) 
}

results <- pairs[, {
  original_file <- Sys.glob(sprintf("data/*/%s.txt", omicspred_id))
  filtered_file <- sprintf("data/omicspred_filtered/%s__%s.txt", protein, endpoint)
  cis_filt_file <- sprintf("data/filtered_cis_pPGS/filtlocus_%s__%s.tsv", protein, endpoint)
  file_cis <- sprintf("data/cis_pPRS/cislocus_%s_%s.locus", protein, omicspred_id,
  trans_filt_file <- sprintf("data/filtered_trans_pPGS/filtlocus_%s__%s.tsv", protein, endpoint))
  
  n_original = count_snps(original_file)
  n_filtered = count_snps(filtered_file)
  n_cis_filtered= count_snps_nohead(cis_filt_file)
  n_cis= count_snps_nohead(file_cis)
  n_trans_filtered= count_snps_nohead(trans_filt_file)
  n_trans_original = n_original - n_cis
  
  #n_cis_original <- get_variants_cis(log_file_cis)
  
  list(
    n_original = n_original,
    n_filtered = n_filtered,
    n_dropped  = n_original - n_filtered,
    n_cis_original = n_cis,
    n_cis_filtered= n_cis_filtered,
    n_cis_dropped = n_cis - n_cis_filtered,
    n_trans_filtered = n_trans_filtered,
    n_trans_original = n_trans_original,
    n_trans_dropped = n_trans_original - n_trans_filtered
  )
}, by = .(protein, endpoint, omicspred_id)]


#write.table(results, "data/results_instrument_selection.tsv",sep="\t", quote=F, row.names = F)

## load
df= fread("data/results_instrument_selection.tsv")

res_filt= df %>% filter(n_dropped>0)
a= merge(res_filt, moresig, by=c("protein","endpoint"))

res_filt[, frac_dropped_cis  := n_cis_dropped / n_cis_original]
res_filt[, frac_dropped_all  := n_dropped / n_original]
res_filt[, frac_dropped_trans  := n_trans_dropped / n_original]

ggplot(res_filt, aes(x = n_dropped, y = n_cis_dropped)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Total SNPs dropped",
    y = "Cis SNPs dropped",
    title = "Are dropped SNPs primarily cis?"
  ) +
  theme_bw()


wilcox.test(
  res_filt$frac_dropped_cis,
  res_filt$frac_dropped_all,
  paired = TRUE
)

### ------ FG filtering 
results <- pairs[, {
  original_file <- Sys.glob(sprintf("data/*/%s.txt", omicspred_id))
  filtered_file <- sprintf("data/FG_omicspred_filtered/%s__%s.txt", protein, endpoint)
  cis_filt_file <- sprintf("data/FG_filtered_cis_pPGS/filtlocus_%s__%s.tsv", protein, endpoint)
  file_cis <- sprintf("data/cis_pPRS/cislocus_%s_%s.locus", protein, omicspred_id,
                      trans_filt_file <- sprintf("data/FG_filtered_trans_pPGS/filtlocus_%s__%s.tsv", protein, endpoint))
  
  n_original = count_snps(original_file)
  n_filtered = count_snps(filtered_file)
  n_cis_filtered= count_snps_nohead(cis_filt_file)
  n_cis= count_snps_nohead(file_cis)
  n_trans_filtered= count_snps_nohead(trans_filt_file)
  n_trans_original = n_original - n_cis
  
  #n_cis_original <- get_variants_cis(log_file_cis)
  
  list(
    n_original = n_original,
    n_filtered = n_filtered,
    n_dropped  = n_original - n_filtered,
    n_cis_original = n_cis,
    n_cis_filtered= n_cis_filtered,
    n_cis_dropped = n_cis - n_cis_filtered,
    n_trans_filtered = n_trans_filtered,
    n_trans_original = n_trans_original,
    n_trans_dropped = n_trans_original - n_trans_filtered
  )
}, by = .(protein, endpoint, omicspred_id)]

write.table(results, "data/FG_results_instrument_selection.tsv",sep="\t", quote=F, row.names = F)

df= fread("data/FG_results_instrument_selection.tsv")

df %>% filter(n_dropped>0 & n_cis_original>0 & n_cis_original==n_cis_dropped)
df %>% filter(n_dropped>0 & n_trans_original>0 & n_trans_original==n_trans_dropped)




