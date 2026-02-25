library(tidyverse)
library(dplyr)
library(data.table)
library(readxl)
library(googlesheets4)
library(ggplot2)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")
# our proteins

ourprot= fread("data/iLR_outliers_prot_pos.tsv")

# their UniProt IDs
olink= fread("data/Olink_trait_validation_results_with_OMICSPRED_ID.csv") %>% select(`OMICSPRED ID`, Gene, `UniProt ID`)
soma= fread("data/Somalogic_trait_validation_results_with_OMICSPRED_ID.csv") %>% select(`OMICSPRED ID`, Gene, `UniProt ID`)

prot= merge(ourprot, olink, by.x = c("protein","omicspred_id"),  by.y=c("Gene","OMICSPRED ID")) # N = 37 
prot2= merge(ourprot, soma, by.x = c("protein","omicspred_id"),  by.y=c("Gene","OMICSPRED ID")) # N = 57

prot= rbind(prot,prot2)

# olink 3K 
# soma 5k or 7K 

# Supplementary Table from https://www.nature.com/articles/s42004-025-01665-1
path = "/users/danifusc/42004_2025_1665_MOESM4_ESM.xlsx"

#sd1 = read_excel(path, sheet = "SD 1 - %CV", skip = 1) %>% as_tibble() %>% left_join(prot, by = c("Analyte" = "UniProt ID"))
sd3 = read_excel(path, sheet = "SD 3- Interplatform Correlation", skip = 1) %>%
  as_tibble() %>% left_join(prot, by =c("UniProt" = "UniProt ID"))

#sd5 = read_excel(path, sheet = "SD 5 - Model (p-value + beta)") %>% as_tibble() %>% left_join(prot_ids, by = "UniProt")
sd8 = read_excel(path, sheet = "SD 8 - Proteins and Platforms") %>% as_tibble() %>% left_join(prot, by =c("UniProt" = "UniProt ID"))
#sd9 = read_excel(path, sheet = "SD 9 - Literature Markers modif") %>% as_tibble() 

# 1 -- How many of our omicspred Olink proteins (n=37) are measured by Olink 5K and SomaScan 
unique(sd3 %>% filter(!is.na(protein) & Platform_1=="Olink 3K") %>% select(protein))
# - Olink 3K = all of them N = 37

unique(sd3 %>% filter(!is.na(protein) & Platform_1=="SomaScan 7K") %>% select(protein))
# - Somalogic 7K = all of them N = 57

unique(sd3 %>% filter(!is.na(protein) & Platform_1=="Olink 5K") %>% select(protein))

unique(sd3 %>% filter(!is.na(protein) & Platform_2=="SomaScan 11K") %>% select(protein))

#  -- Are the protein measurements correlated? 

### proteins present in all 3 platforms Olink 3K, SomaScan 7K and 11K = 93
common_prot= unique(sd3 %>% filter(!is.na(protein) & Platform_2=="SomaScan 11K") %>% select(protein))

### restricted to our 93 proteins 
bind_rows(
  sd3 %>% filter(!is.na(protein)) %>% mutate(`Spearman Rho` = as.numeric(`Spearman Rho`)),
  sd3 %>% filter(!is.na(protein)) %>% mutate(`Spearman Rho` = as.numeric(`Spearman Rho`)) %>%
    mutate(temp = Platform_1, Platform_1 = Platform_2, Platform_2 = temp) %>% select(-temp)
) %>% 
  filter(Platform_1 == "Olink 3K") %>% 
  filter(Platform_2 %in% c("SomaScan 11K", "SomaScan 7K")) %>%
  filter(protein %in% common_prot$protein) %>%
  group_by(Platform_2) %>% 
  summarize(median_rho = median(`Spearman Rho`), mean_rho = mean(`Spearman Rho`))

### all proteins available in Olink 3K, SomaScan 11K and SomaScan 7K
sd3_full = read_excel(path, sheet = "SD 3- Interplatform Correlation", skip = 1) %>%
  as_tibble() %>% filter(Platform_1 == "Olink 3K") %>% 
  filter(Platform_2 %in% c("SomaScan 11K", "SomaScan 7K")) %>%
  group_by(Platform_2_Analyte) %>%
  filter(n_distinct(Platform_2) == 2) %>%
  ungroup()

## median Rho
bind_rows(
sd3_full %>% mutate(`Spearman Rho` = as.numeric(`Spearman Rho`)),
sd3_full %>% mutate(`Spearman Rho` = as.numeric(`Spearman Rho`)) %>%
  mutate(temp = Platform_1, Platform_1 = Platform_2, Platform_2 = temp) %>% select(-temp)
) %>% 
  filter(Platform_1 == "Olink 3K") %>% 
  filter(Platform_2 %in% c("SomaScan 11K", "SomaScan 7K")) %>%
  group_by(Platform_2) %>% 
  summarize(median_rho = median(`Spearman Rho`), mean_rho = mean(`Spearman Rho`))


comp= bind_rows(
  sd3 %>% filter(!is.na(protein)) %>% mutate(`Spearman Rho` = as.numeric(`Spearman Rho`)),
  sd3 %>% filter(!is.na(protein)) %>% mutate(`Spearman Rho` = as.numeric(`Spearman Rho`)) %>% mutate(temp = Platform_1, Platform_1 = Platform_2, Platform_2 = temp)
) %>% 
  filter(Platform_1 == "Olink 3K") %>% 
  filter(Platform_2 %in% c("Olink 5K", "SomaScan 7K", "SomaScan 11K")) %>% 
  group_by(protein, Platform_2) %>% 
  summarize(`Spearman Rho` = median(`Spearman Rho`, na.rm = FALSE)) %>%
  ungroup() %>% 
  complete(protein, Platform_2, fill = list(`Spearman Rho` = NA)) %>% 
  
  ggplot(aes(y = Platform_2, x = protein, fill = `Spearman Rho`)) + 
  geom_tile() + 
  labs(x = "Spearman Rho with Olink 3K", y = "", fill = "Spearman correlation between proteins") + 
  scale_fill_gradient2(low = "#1A85FF", high = "#D41159", mid = "white", midpoint = 0,
                       na.value = "lightgray", limits = c(-1, 1)) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

pdf("plots/protein_comparison_platforms.pdf", width = 15.5, height = 4.5)
comp
dev.off()

png("plots/protein_comparison_platforms.png", width = 15.5, height = 4.5, units = "in", res = 300)
comp
dev.off()


# What is the distribution of our proteins correlations (olink 3K - SomaScan 7K) vs all their proteins ? 

their_prot= sd3_full %>% filter(Platform_2=="SomaScan 7K") %>%
  mutate(Source="All proteins available in\n Kirsher et al.")

our_prot= sd3 %>% filter(Platform_1 == "Olink 3K" & Platform_2=="SomaScan 7K" & protein %in% common_prot$protein)  %>%
  mutate(Source="Restricted to our proteins\n (N=93)")

df_long <- bind_rows(their_prot, our_prot)


distplot= ggplot(df_long, aes(x = as.numeric(`Spearman Rho`),
                    group = Source,
                    fill=Source)) +
  scale_fill_manual(values= c("All proteins available in\n Kirsher et al."="#4DA1D0",
                               "Restricted to our proteins\n (N=93)"="#E69F00")) +
  geom_density(alpha = 0.6) +
  labs(x = "Spearman Rho between Olink 3K and SomaScan 7K", y = "Density") +
  theme_classic() +
  theme(legend.title = element_blank())

pdf("plots/distribution_comparison_platforms.pdf", width = 11, height = 8)
distplot
dev.off()

png("plots/distribution_comparison_platforms.png", width = 11, height = 8, units = "in", res = 300)
distplot
dev.off()

