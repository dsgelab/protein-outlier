library(biomaRt)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

library(dplyr)
library(readr)


dir= "/scratch/project_2007428/projects/prj_100_pprs_discordance/data/"
prot= fread("data/idefix_olinksoma/pheno_residuals.tsv")
cis_prs= fread("data/cis_genoresiduals.tsv")

### compute cis-R2 of the proteins
r2= function(dfprs) {
  common_prot = intersect(names(prot[,-1]), names(dfprs[,-1]))
  results=  map(common_prot, function(col) {
    model= lm(prot[[col]] ~ dfprs[[col]])
    summ = summary(model)
    tibble(
      protein = col,
      beta = coef(model)[2],
      se = coef(summ)[2, "Std. Error"],
      pval = coef(summ)[2, "Pr(>|t|)"],
      r2 = summ$r.squared,
      adj_r2 = summ$adj.r.squared
    )
  })
  results_df = bind_rows(results)
}

cis_r2= r2(cis_prs)


cis_window_bp= 1e06

#### cis loci
loci = fread("data/iLR_outliers_prot_pos.tsv") %>%
  dplyr::select(!R2) %>% 
  as_tibble() %>%
  mutate(
    folder = case_when(
      type %in% c("olink", "Olink") ~ "Olink",
      type %in% c("somalogic", "SomaScan", "somascan") ~ "SomaScan",
      TRUE ~ NA_character_
    ),
    file = file.path("data", folder, paste0(omicspred_id, ".txt")),
    cis_start = pmax(1L, as.integer(start_position - cis_window_bp)),
    cis_end   = as.integer(end_position + cis_window_bp),
    chromosome_name = as.character(chromosome_name)
  )

process_one = function(row) {
  # row is a 1-row tibble
  f <- row$file[[1]]
  if (!file.exists(f)) {
    return(tibble(
      protein = row$protein, omicspred_id = row$omicspred_id, type = row$type, R2 = row$R2,
      top_rsid = NA_character_, top_effect_weight = NA_real_, top_chr = NA_character_, top_pos = NA_integer_
    ))
  }

  
  w <- fread(f) %>% as_tibble()
  
  
  # columns  rsid, chr_name, chr_position, effect_weight
  w_cis <- w %>%
    mutate(
      chr = gsub("^chr", "", as.character(chr_name)),
      pos = as.integer(chr_position)
    ) %>%
    filter(
      chr == as.character(row$chromosome_name[[1]]),
      pos >= row$cis_start[[1]],
      pos <= row$cis_end[[1]]
    )
  
  if (nrow(w_cis) == 0) {
    return(tibble(
      protein = row$protein, omicspred_id = row$omicspred_id, type = row$type,
      top_rsid = NA_character_, top_effect_weight = NA_real_, top_chr = NA_character_, top_pos = NA_integer_
    ))
  }
  
  w_cis %>%
    slice_max(order_by = abs(effect_weight), n = 1, with_ties = FALSE) %>%
    transmute(
      protein = row$protein,
      omicspred_id = row$omicspred_id,
      type = row$type,
      top_rsid = rsid,
      top_effect_weight = effect_weight,
      top_chr = chr,
      top_pos = pos
    )
}

top_cis <- loci %>%
  split(.$omicspred_id) %>%
  map_dfr(~process_one(.x[1, ]))
  


mart <- biomaRt::useEnsembl(
  biomart = "snp",
  dataset = "hsapiens_snp",
  mirror  = "useast",
  GRCh = 37 
)

get_conseq <- function(x) {
  biomaRt::getBM(
    attributes = c("refsnp_id", "consequence_type_tv"),
    filters    = "snp_filter",   # if this errors, switch to filters="refsnp_id"
    values     = x,
    mart       = mart
  )
}

batch_size <- 500
batches <- split(top_cis$top_rsid, ceiling(seq_along(top_cis$top_rsid) / batch_size))

anno_raw <- map_dfr(batches, ~{
  tryCatch(get_conseq(.x), error = function(e) {
    message("Batch failed: ", conditionMessage(e))
    tibble(refsnp_id = character(), consequence_type_tv = character())
  })
})


# 4) flag PAV-like consequences
pav_terms <- c(
  "missense_variant",
  "stop_gained","stop_lost","start_lost",
  "frameshift_variant",
  "inframe_insertion","inframe_deletion",
  "splice_donor_variant","splice_acceptor_variant"
)
pav_regex <- paste0("\\b(", paste(pav_terms, collapse="|"), ")\\b")

anno_flag <- anno_raw %>%
  mutate(pav_hit = str_detect(consequence_type_tv, pav_regex)) %>%
  group_by(refsnp_id) %>%
  summarise(
    PAV_like = any(pav_hit),
    consequences = paste(sort(unique(consequence_type_tv)), collapse = ","),
    .groups = "drop"
  ) %>%
  rename(top_rsid = refsnp_id)

top_cis_annot <- top_cis %>%
  left_join(anno_flag, by = "top_rsid") %>%
  mutate(
    PAV_like = if_else(is.na(PAV_like), NA, PAV_like),
    consequences = if_else(is.na(consequences), "", consequences)
  )


top_cis_annot %>%
  count(PAV_like)

top_cis_annot = merge(top_cis_annot, cis_r2, by= "protein")

m1 <- glm(PAV_like ~ r2 + type, data=top_cis_annot, family=binomial())
summary(m1)

# Tail enrichment: top 5% R2
cut <- quantile(top_cis_annot$r2, 0.95, na.rm=TRUE)
df <- top_cis_annot %>% mutate(R2_top5 = as.integer(r2 >= cut))
m2 <- glm(PAV_like ~ R2_top5 + type, data=df, family=binomial())
summary(m2)


### see how many proteins are present in the paper Suhre et al. 
### How many of them are PAV-like : 
sd3= fread("data/SD3_paperSuhre_Epitopes.txt", skip = 2)
our= fread("data/R2_94proteins_olinksoma.tsv")
### common proteins 
common= intersect(sd3$SYMBOL, our$protname)

### how many of them have been tested by Suhre et al.
sd2=  readxl::read_xlsx("~/SD2_Suhre.xlsx", skip = 2)
common_all= intersect(sd2$SYMBOL, our$protname )

