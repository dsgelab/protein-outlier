library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(survival)
library(ggplot2)

setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

# ── Load data ──────────────────────────────────────────────
allproteins   <- fread("data/idefix_olinksoma/pheno_residuals.tsv")
setDT(allproteins[, (names(allproteins)[-1]) := lapply(.SD, scale), .SDcols = names(allproteins)[-1]])

resid_nocov   <- fread("data/idefix_olinksoma/residuals_output_nocov.tsv")
resid_nocov$AggAbsResid <- NULL
setDT(resid_nocov[, (names(resid_nocov)[-1]) := lapply(.SD, scale), .SDcols = names(resid_nocov)[-1]])
colnames(resid_nocov) <- gsub("StandResid_", "", colnames(resid_nocov))


intervene_fi  <- arrow::read_parquet("/scratch/project_2007428/projects/prj_010_phenotype_file/share/ukb78537_INTERVENE_phenotype-2024-08-01.parquet")
intervene_fi  <- intervene_fi %>% dplyr::select(ID, SEX, END_OF_FOLLOWUP, contains("_DATE"))
intervene_fi  <- intervene_fi %>% filter(ID %in% allproteins$eid)

phenos        <- fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv")
selphenos     <- phenos %>% dplyr::select(ID = eid, age = '21022-0.0', sex = '31-0.0',
                                   baseline = '53-0.0', death = '40000-0.0') %>%
  filter(ID %in% intervene_fi$ID)

intervene_fi  <- merge(intervene_fi, selphenos, by = "ID")

pairs_to_test= fread("data/results_logreg_all.tsv")
pairs_to_test= pairs_to_test %>%  mutate(p_adjprot= p.adjust(pval_prot, method="fdr"),
                          p_adjresid= p.adjust(pval_resid, method="fdr")) %>% 
  filter((p_adjprot < 0.05 | p_adjresid < 0.05)) %>% dplyr::select(protein, endpoint)


# ── Original incident_logif function ──────────────────────────────
incident_logif <- function(prot_df, idcol, intervene_data) {
  allendpoints <- data.frame()
  for (dis in colnames(intervene_data %>% dplyr::select(contains("_DATE")))) {
    
    df <- intervene_data %>% dplyr::select(ID, END_OF_FOLLOWUP, dis, baseline, age, sex)
    
    # compute n_cases once per endpoint
    df_bin_check <- df %>%
      mutate(fol = difftime(!!as.name(dis), baseline, units = "weeks") / 52.25) %>%
      filter(fol > 0.5 | is.na(fol)) %>%
      mutate(!!as.name(dis) := ifelse(is.na(!!as.name(dis)), 0, 1))
    
    n_cases_endpoint <- sum(df_bin_check[[dis]] == 1, na.rm = TRUE)
    
    message("dis: ", dis, " | n_cases: ", n_cases_endpoint, " | n_total: ", nrow(df_bin_check))
    
    # ── Skip endpoint if fewer than 10 cases ─────────────────────────────
    if (n_cases_endpoint < 10) {
      message("SKIPPING — fewer than 10 cases for dis: ", dis)
      next
    }
    
    allresults <- data.frame()
    
    proteins_to_run <- pairs_to_test %>%
      filter(endpoint == gsub("_DATE","",dis)) %>%
      pull(protein)
    
    if (length(proteins_to_run) == 0) next
    
    for (prot in proteins_to_run) {
      p        <- prot_df %>% dplyr::select(ID = idcol, prot)
      df_merge <- merge(df, p, by = "ID")
      df_merge <- df_merge %>%
        mutate(fol = difftime(!!as.name(dis), baseline, units = "weeks") / 52.25)
      df_filt  <- df_merge %>% filter(fol > 0.5 | is.na(fol))
      df_bin   <- df_filt %>%
        mutate(!!as.name(dis) := ifelse(is.na(!!as.name(dis)), 0, 1))
      
      n_cases <- sum(df_bin[[dis]] == 1, na.rm = TRUE)
      n_total <- nrow(df_bin)
      
      formula   <- as.formula(paste0(dis, "~", prot, " + age + sex"))
      
      glm_model <- tryCatch({
        glm(formula, data = df_bin, family = "binomial")
      }, error = function(e) {
        message("ERROR — dis: ", dis, " | prot: ", prot,
                " | n_cases: ", n_cases, " | msg: ", e$message)
        return(NULL)
      }, warning = function(w) {
        message("WARNING — dis: ", dis, " | prot: ", prot,
                " | n_cases: ", n_cases, " | msg: ", w$message)
        return(suppressWarnings(glm(formula, data = df_bin, family = "binomial")))
      })
      
      if (is.null(glm_model)) next
      
      glmdf         <- as.data.frame(summary(glm_model)$coefficients[2, , drop = FALSE])
      glmdf$protein <- prot
      glmdf$n_cases <- n_cases
      glmdf$n_total <- n_total
      
      allresults <- rbind(glmdf, allresults)
    }
    
    if (nrow(allresults) == 0) {
      message("NO RESULTS for dis: ", dis)
      next
    }
    
    allresults$endpoint <- dis
    allendpoints        <- rbind(allresults, allendpoints)
  }
  return(allendpoints)
}
# ── Down-sampling setup ─────────────────────────────────────
all_ids <- intervene_fi$ID ## ids 
full_N <- length(intervene_fi$ID) ## full sample size

# Define sample sizes

sample_sizes <- c(1000, 2000, 5000, 10000, 20000, 30000, full_N)
sample_sizes <- sort(sample_sizes[sample_sizes <= full_N])


# ── Down-sampling loop ─────────────────────────────────────────────────────────
all_downsampling_results <- list()

for (N in sample_sizes) {
  cat("Running N =", N, "\n")
  
  set.seed(1234)
  sub_ids <- if (N == full_N) all_ids else sample(all_ids, N, replace = FALSE)
  
  sub_proteins  <- allproteins  %>% filter(eid   %in% sub_ids)
  sub_resid     <- resid_nocov  %>% filter(IndID %in% sub_ids)
  sub_intervene <- intervene_fi %>% filter(ID    %in% sub_ids)
  
  res_prot <- incident_logif(sub_proteins, "eid", sub_intervene) %>%
    mutate(model = "raw_protein", N = N)
  
  res_resid <- incident_logif(sub_resid, "IndID", sub_intervene) %>%
    mutate(model = "ppgs_adjusted", N = N)
  
  all_downsampling_results[[as.character(N)]] <- bind_rows(res_prot, res_resid)
  
}

# ── Combine all results ────────────────────────────────────────────────────────
ds_results <- bind_rows(all_downsampling_results)
ds_results$protein = gsub("StandResid_","",ds_results$protein)

# Rename columns
ds_results <- ds_results %>%
  rename(estimate = Estimate, SE = `Std. Error`, zvalue = `z value`, pval = `Pr(>|z|)`)

write.table(ds_results, file = "data/downsampling_all_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── Summarise stability metrics across iterations ──────────────────────────────
ds_results= fread("data/downsampling_all_results.tsv")

# Find pairs present at ALL sample sizes
ds_results_wide <- ds_results %>%
  pivot_wider(
    names_from  = model,
    values_from = c(estimate, SE, zvalue, pval)
  ) %>%
  select(-n_cases, -n_total) %>%  # same values, drop duplicate
  mutate(
    p_adj_raw  = p.adjust(pval_raw_protein,      method = "fdr"),
    p_adj_ppgs = p.adjust(pval_ppgs_adjusted, method = "fdr"),
    sig_either = p_adj_raw < 0.05 | p_adj_ppgs < 0.05
  )

ds_results_wide_sig <- ds_results_wide %>%
  group_by(protein, endpoint) %>%
  filter(
    n_distinct(N) == n_distinct(ds_results_wide$N),  # present at all N
    all(sig_either, na.rm = TRUE)                     # significant at all N
  ) %>%
  ungroup()

# Sanity check
ds_results_wide_sig %>%
  group_by(N) %>%
  summarise(n_pairs = n()) %>%
  print()


plot= ds_results_wide_sig %>%
  mutate(N_label = factor(N,
                          levels = c(1000, 2000, 5000, 10000, 20000, 30000, 39871),
                          labels = c("1K", "2K", "5K", "10K", "20K", "30K", "Full"))) %>%
  ggplot(aes(x = zvalue_raw_protein, y = zvalue_ppgs_adjusted)) +
  geom_point(shape = 21, size = 1.5, fill = "lightgrey", colour = "black", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey56") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey56") +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.3,
              se = FALSE, inherit.aes = FALSE,
              aes(x = zvalue_raw_protein, y = zvalue_ppgs_adjusted)) +
  facet_wrap(~ N_label, nrow = 2) +
  labs(x = "Association of unadjusted protein with incident disease (z-score)",
       y = "Association of genetically adjusted protein with incident disease (z-score)") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9),
        strip.text   = element_text(face = "bold", size = 11))

ggsave(
  "plots/Supplementary_downsampling.pdf",
  plot = plot,
  device = cairo_pdf,
  width = 180,
  height = 150,   # adapt to figure
  units = "mm"
)

ggsave(
  "plots/Supplementary_downsampling.png",
  plot = plot,
  width = 180,
  height = 150,
  units = "mm",
  dpi = 300
)







