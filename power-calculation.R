# Power calculation 
# Ratio of the required sample sizes R = (z_prot / z_resid)^2

library("tidyverse")
library("data.table")

df = fread("data/results_logreg_all.tsv") %>% as_tibble() 

# Summary for FDR-significant protein-disease pairs
df %>% 
  filter(p.adjust(pval_prot, "fdr") < 0.05) %>% 
  mutate(R = (zvalue_prot / zvalue_resid)^2) %>% 
  pull(R) %>% 
  summary()

# Supplementary table 
df %>% 
  filter(p.adjust(pval_prot, "fdr") < 0.05) %>% 
  mutate(R = (zvalue_prot / zvalue_resid)^2) %>% 
  select(protein, endpoint, estimate_prot, SE_prot, estimate_resid, SE_resid, R) %>% 
  arrange(R) %>% 
  fwrite("results/power-calculation.csv")
