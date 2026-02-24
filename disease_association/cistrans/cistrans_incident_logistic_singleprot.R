library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(survival)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(purrr)
library(ggrepel)
setwd("/scratch/project_2007428/projects/prj_100_pprs_discordance/")

prot= fread("data/idefix_olinksoma/pheno_residuals.tsv")
#scale
setDT(prot[, (names(prot)[-1]) := lapply(.SD, scale), .SDcols = names(prot)[-1]])

cis_resid= fread("data/cis_residuals_output_nocov.tsv") %>% select(!AggAbsResid)
setDT(cis_resid[, (names(cis_resid)[-1]) := lapply(.SD, scale), .SDcols = names(cis_resid)[-1]])

trans_resid= fread("data/trans_residuals_output_nocov.tsv") %>% select(!AggAbsResid)
setDT(trans_resid[, (names(trans_resid)[-1]) := lapply(.SD, scale), .SDcols = names(trans_resid)[-1]])

cis_prs= fread("data/cis_genoresiduals.tsv")
setDT(cis_prs[, (names(cis_prs)[-1]) := lapply(.SD, scale), .SDcols = names(cis_prs)[-1]])

trans_prs= fread("data/trans_genoresiduals.tsv")
setDT(trans_prs[, (names(trans_prs)[-1]) := lapply(.SD, scale), .SDcols = names(trans_prs)[-1]])

### file with first occurrences of INTERVENE phenotypes 
intervene_fi= arrow::read_parquet("/scratch/project_2007428/projects/prj_010_phenotype_file/share/ukb78537_INTERVENE_phenotype-2024-08-01.parquet")
intervene_fi= intervene_fi %>% select(ID,SEX,END_OF_FOLLOWUP, contains("_DATE"))
### filter for our individuals
intervene_fi = intervene_fi %>% filter(ID %in% prot$eid)

### filter phenotypes for the intervene. note that from the 39 endpoints, BMI and Covid-19 severity are not present as first occurrences, so the final number is 37 
definitions= fread("data/Intervene_flagship_endpoint_collection_Definitions.tsv") 
definitions = definitions %>% select(Endpoint,`FinnGen endpoint`) 

## load age and sex covariates, first assessment visit date and death
phenos=fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv")

selphenos= phenos %>% select(ID=eid, age='21022-0.0',
                             sex='31-0.0', baseline='53-0.0', death='40000-0.0')
selphenos = selphenos %>% filter(ID %in% intervene_fi$ID)

#### The column "END_OF_FOLLOWUP" is the age of the first record of a disease diagnosis (for individuals with the diseases), age at death for a cause other than the disease, age at last record available in the registries or electronic health records or age 80, whatever happened first
intervene_fi= merge(intervene_fi, selphenos, by="ID")

################################# INCIDENT LOGISTIC REGRESSION
incident_logif= function(prot_df, idcol) {
  allendpoints = data.frame()
  for (dis in colnames(intervene_fi %>% select(contains("_DATE")))) {
    df= intervene_fi %>% select(ID, END_OF_FOLLOWUP, dis, baseline, age, sex) 
    allresults= data.frame()
    for (prot in colnames(prot_df)[-1]) {
      p= prot_df %>% select(ID= idcol, prot)
      df_merge= merge(df, p, by="ID")
      df_merge= df_merge %>% mutate(fol = difftime(!!as.name(dis), baseline, units = "weeks")/52.25) ###note that I converted in years!!
      
      ## excluding all prevalent cases or incident cases recorded within the first 6 months of follow-up
      df_filt= df_merge %>% filter(fol>0.5 | is.na(fol))
      
      ### binary outcome, 1= incidence, 0 = right censored
      df_bin= df_filt %>% mutate(!!as.name(dis) := ifelse(is.na(!!as.name(dis)),0,1))
      
      ### Logistic regression
      formula= as.formula(paste0(dis,"~",prot," + age + sex"))
      glm_model = glm(formula, data = df_bin, family = "binomial")
      
      glmdf= as.data.frame(summary(glm_model)$coefficients[2, , drop = FALSE])
      glmdf$protein= prot
      
      allresults = rbind(glmdf, allresults)
      
    }
    allresults$endpoint= dis
    allendpoints= rbind(allresults, allendpoints)
  }
  return(allendpoints)
}

cis_lr= incident_logif(cis_resid, "IndID")
cis_lr$protein = gsub("StandResid_","", cis_lr$protein)
cis_lr= cis_lr %>% mutate(endpoint= gsub("_DATE","",endpoint))
cis_lr$Model= "cis-residuals"

cisprs_lr= incident_logif(cis_prs, "eid")
cisprs_lr$Model= "cis-PRS"

transprs_lr= incident_logif(trans_prs, "eid")
transprs_lr$Model= "trans-PRS"

trans_lr= incident_logif(trans_resid, "IndID")
trans_lr$protein = gsub("StandResid_","", trans_lr$protein)
trans_lr= trans_lr %>% mutate(endpoint= gsub("_DATE","",endpoint))
trans_lr$Model= "trans-residuals"

inc_resid= fread("data/incident_logreg_singleresidualsnocov.tsv") %>% mutate(endpoint= gsub("_DATE","",endpoint)) %>%
  rename(Estimate=estimate_resid, `Std. Error`=SE_resid,
                                  `z value`=zvalue_resid, `Pr(>|z|)`=pval_resid)
inc_resid$Model= "residuals"

inc_prot= fread("data/incident_logreg_singleproteins.tsv") %>% mutate(endpoint= gsub("_DATE","",endpoint)) %>%
  rename(Estimate=estimate_prot, `Std. Error`=SE_prot,
         `z value`=zvalue_prot, `Pr(>|z|)`=pval_prot)
inc_prot$Model= "protein"

inc_prs= fread("data/incident_logreg_singlePRS.tsv") %>% mutate(endpoint= gsub("_DATE","",endpoint)) %>%
  rename(Estimate=estimate_PRS, `Std. Error`=SE_PRS,
         `z value`=zvalue_PRS, `Pr(>|z|)`=pval_PRS)
inc_prs$Model= "PRS"


all_lr= rbind(cis_lr,trans_lr,inc_resid,inc_prot, inc_prs, transprs_lr, cisprs_lr)
all_lr$endpoint= gsub("_DATE","",all_lr$endpoint)
all_lr= merge(all_lr,definitions, by.x="endpoint",by.y="FinnGen endpoint")

all_lr$Model= factor(all_lr$Model, levels = c("protein", "PRS","residuals", "cis-PRS","trans-PRS","cis-residuals", "trans-residuals"))
### ODDS RATIO  
all_lr = all_lr %>% mutate(OR= exp(Estimate),
                           lower_ci= exp(Estimate - 1.96*`Std. Error`),
                           upper_ci=exp(Estimate + 1.96*`Std. Error`))
write.table(all_lr, file="data/cistransresiduals_logisticregressions.tsv", row.names = F,sep="\t", quote=F)

########################################### LOAD 
cis_resid= fread("data/cis_residuals_output_nocov.tsv") %>% select(!AggAbsResid) %>% rename(eid=IndID)
colnames(cis_resid)= gsub("StandResid_","",colnames(cis_resid))
trans_resid= fread("data/trans_residuals_output_nocov.tsv") %>% select(!AggAbsResid) %>% rename(eid=IndID)
colnames(trans_resid)= gsub("StandResid_","",colnames(trans_resid))
resid= fread("data/idefix_olinksoma/residuals_output.tsv") %>% select(!AggAbsResid) %>% rename(eid=IndID)
colnames(resid)= gsub("StandResid_","",colnames(resid))

cis_prs= fread("data/cis_genoresiduals.tsv")
trans_prs= fread("data/trans_genoresiduals.tsv")
prs= fread("data/idefix_olinksoma/geno_residuals.tsv") 


##############################
dir= "/scratch/project_2007428/projects/prj_100_pprs_discordance/data/"
prot= fread("data/idefix_olinksoma/pheno_residuals.tsv")
cis_prs= fread("data/cis_genoresiduals.tsv")
trans_prs= fread("data/trans_genoresiduals.tsv")

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
trans_r2= r2(trans_prs)

all_r2= merge(cis_r2,trans_r2, by="protein",suffixes = c("_cis","_trans"))

r2p= ggplot(data= all_r2, aes(x= r2_cis, y= r2_trans)) +
  geom_point(size = 2.5) +
  theme_minimal(base_size=13) + 
  #geom_smooth(method = "lm", colour = "grey56", linewidth = 0.4, se = FALSE, alpha=0.5) +
  geom_abline(slope=1, intercept=0, color = "black", linetype ="dashed") +
  labs(x= "R2 (protein ~ cis-PRS)", y="R2 (protein ~ trans-PRS)") +
  xlim(0,0.82) + ylim(0,0.82) + 
  geom_text_repel(aes(label = protein),
                  size=2.5 , vjust = -1, hjust=-0.5, max.overlaps = 20)

#### LOAD INCIDENT LOGISTIC REGRESSIONS
inc_prot= fread("data/incident_logreg_singleproteins.tsv")
inc_cistrans= fread("data/cistransresiduals_logisticregressions.tsv")
inc_cistrans = inc_cistrans %>% rename(zscore = `z value`,
                                       pval= "Pr(>|z|)",
                                       se= "Std. Error")

cistrans_wide = inc_cistrans %>%
  pivot_wider(
    id_cols = c(endpoint, protein, Endpoint),  # keep these as identifiers
    names_from = Model,
    values_from = c(Estimate, se, zscore, pval, OR, lower_ci, upper_ci),
    names_glue = "{.value}_{Model}"
  )

cistrans_wide= cistrans_wide %>% mutate(p_adjprot= p.adjust(pval_protein, method="fdr"),
                                        p_adjresid= p.adjust(pval_residuals, method="fdr"))


#### beta difference ( residuals - protein )
cis_wide = cistrans_wide %>% mutate(
  `cis-diffresprot` = `Estimate_cis-residuals`-Estimate_protein,
  `cis-se_diffresprot`= sqrt(`se_cis-residuals`^2 + se_protein^2),
  `cis-z_diffresprot` = `cis-diffresprot`/`cis-se_diffresprot`,
  `cis-p_diffresprot` = 2 * (1 - pnorm(abs(`cis-z_diffresprot`)))
)
cis_wide= merge(cis_wide, cis_r2, by="protein")

cis_wide = cis_wide %>% mutate(p_adjcisres= p.adjust(`pval_cis-residuals`, method = "fdr"),
                               sig_fdr= ifelse(p_adjcisres<0.05 | p_adjprot<0.05, "yes","no"),
                               sig_zdiff= ifelse(`cis-p_diffresprot`<0.05, "yes","no"))
#cis_wide_sigdiff= cis_wide %>% filter(`cis-p_diffresprot`<0.05)

trans_wide = cistrans_wide %>% mutate(
  `trans-diffresprot`= `Estimate_trans-residuals`-Estimate_protein,
  `trans-se_diffresprot`= sqrt(`se_trans-residuals`^2 + se_protein^2),
  `trans-z_diffresprot` = `trans-diffresprot`/`trans-se_diffresprot`,
  `trans-p_diffresprot`= 2 * (1 - pnorm(abs(`trans-z_diffresprot`)))
)
trans_wide= merge(trans_wide, trans_r2, by="protein")
trans_wide = trans_wide %>% mutate(p_adjtransres= p.adjust(`pval_trans-residuals`, method = "fdr"),
                                   sig_fdr= ifelse(p_adjtransres<0.05 | p_adjprot<0.05, "yes","no"),
                                   sig_zdiff= ifelse(`trans-p_diffresprot`<0.05, "yes","no"))

cis_wide_padj = cis_wide %>% filter(p_adjprot<0.05 | p_adjcisres<0.05)
trans_wide_padj = trans_wide %>% filter(p_adjprot<0.05 | p_adjtransres<0.05)

##### SUPPLEMENTARY FIGURE DIFFERENCE IN Z-SCORES VS R2 
p1= ggplot(data= cis_wide_padj,
           aes(x= r2, y= abs(`cis-z_diffresprot`), 
               fill= ifelse(`cis-p_diffresprot`<0.05,"yes","no"))
) +
  geom_point(shape=21, size = 2.5, alpha=0.8) +
  scale_fill_manual(values = c("no" = "lightgrey", "yes" = "blue4"),
                    name = "Significant difference between\ncis-adjusted and unadjusted\nprotein-disease pairs (P < 0.05)") +
  geom_point(shape=21,size = 2.5, alpha=0.8) +
  theme_minimal() + 
  xlim(0,0.85) +
  ylim(0,9) +
  theme(
    axis.title.x = element_text(face = "bold",size = 11),
    axis.title.y = element_text(face = "bold",size = 10),
    plot.tag = element_text(face = "bold", size = 14)) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.4, se = FALSE,
              inherit.aes = FALSE,
              aes(x = r2, y = abs(`cis-z_diffresprot`))) +
  labs(x= bquote(bold(R^2 ~ "(unadjusted protein ~ cis-PGS)")),
       y="Difference between cis-adjusted\n and unadjusted protein ( |z-score| )")

p2= ggplot(data= trans_wide_padj,
           aes(x= r2, y= abs(`trans-z_diffresprot`),
               fill= ifelse(`trans-p_diffresprot`<0.05,"yes","no"))
) +
  geom_point(shape=21, size = 2.5, alpha=0.8) +
  scale_fill_manual(values = c("no" = "lightgrey", "yes" = "blue4"),
                    name = "Significant difference between\ntrans-adjusted and unadjusted\nprotein-disease pairs (P < 0.05)") +
  theme_minimal() + 
  xlim(0,0.85) +
  ylim(0,9) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.4, se = FALSE,  
              inherit.aes = FALSE,
              aes(x = r2, y = abs(`trans-z_diffresprot`))) +
  theme(
    axis.title.x = element_text(face = "bold",size = 11),
    axis.title.y = element_text(face = "bold",size = 10),
    plot.tag = element_text(face = "bold", size = 14)) +
  labs(x= bquote(bold(R^2 ~ "(unadjusted protein ~ trans-PGS)")),
       y="Difference between trans-adjusted\n and unadjusted protein ( |z-score| )")

pdf("plots/Supplementary_FigDIFFR2.pdf", height=9.2, width = 7.8)
(p1 / p2) + plot_annotation(tag_levels = "A")
dev.off()

png("plots/Supplementary_FigDIFFR2.png", height=9.2, width = 7.8, units = "in", res = 300)
(p1 / p2) + plot_annotation(tag_levels = "A")
dev.off()



#####################
cis_wide_padj= cis_wide %>% filter(p_adjcisres<0.05 | p_adjprot<0.05)
nrow(cis_wide_padj %>% filter(`cis-p_diffresprot`<0.05 & abs(`zscore_cis-residuals`)>abs(zscore_protein)))
nrow(cis_wide_padj %>% filter(`cis-p_diffresprot`<0.05 & abs(`zscore_cis-residuals`)<abs(zscore_protein)))


trans_wide_padj= trans_wide %>% filter(p_adjtransres<0.05 | p_adjprot<0.05)
nrow(trans_wide_padj %>% filter(`trans-p_diffresprot`<0.05 & abs(`zscore_trans-residuals`)>abs(zscore_protein)))
nrow(trans_wide_padj %>% filter(`trans-p_diffresprot`<0.05 & abs(`zscore_trans-residuals`)<abs(zscore_protein)))

############# SUPPLEMENTARY FIGURE LOGISTIC REGRESSIONS CIS TRANS
p <- function(data, x, y, sig_fdr, fill, label, legendpos="right") {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_point(aes(shape = .data[[sig_fdr]], fill = .data[[fill]]), size = 2) +
    
    scale_fill_manual(values = c("no" = "lightgrey", "yes" = "blue4"), 
                      name = "Significant difference between\ntrans-adjusted and unadjusted\nprotein-disease pairs (P < 0.05)") + 
    scale_shape_manual(values = c("yes" = 21, "no" = 24), 
                       name = "FDR < 0.05 in at least\none model") +
    
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey56") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey56") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey56") +
    
    geom_smooth(data = data, aes(x = .data[[x]], y = .data[[y]]), 
                method = "lm", colour = "black", linewidth = 0.3, se = FALSE, inherit.aes = FALSE) +
    
    labs(
      x = "Association of protein with incident disease (z-score)",
      y = paste0("Association of ", label, "\nwith incident disease (z-score)")
    ) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(face = "bold", size = 11),
      axis.title.y = element_text(face = "bold", size = 11),
      legend.position = legendpos,
      plot.tag = element_text(face = "bold", size = 14)
    )
}

#### create plot axes limits so they are the same 
rng_x = range(c(cis_wide$zscore_protein, trans_wide$zscore_protein))
rng_y = range(c(cis_wide$`zscore_cis-residuals`, trans_wide$`zscore_trans-residuals`))

  
pdf("plots/SupplementaryFigLRcistrans.pdf", width = 15, height = 7)
(p(cis_wide, "zscore_protein", "zscore_cis-residuals", "sig_fdr", "sig_zdiff", "cis-genetically-adjusted protein ", "none") +
p(trans_wide, "zscore_protein", "zscore_trans-residuals", "sig_fdr", "sig_zdiff", "trans-genetically-adjusted protein ") +
  plot_annotation(tag_levels = "A")) & coord_cartesian(xlim = rng_x, ylim = rng_y)
dev.off()

png("plots/SupplementaryFigLRcistrans.png", width = 15, height = 7, units = "in", res = 300)
(p(cis_wide, "zscore_protein", "zscore_cis-residuals", "sig_fdr", "sig_zdiff", "cis-genetically-adjusted protein ", "none") +
    p(trans_wide, "zscore_protein", "zscore_trans-residuals", "sig_fdr", "sig_zdiff", "trans-genetically-adjusted protein ") +
    plot_annotation(tag_levels = "A")) & coord_cartesian(xlim = rng_x, ylim = rng_y)
dev.off()



##########################################
outlierstrans= trans_wide %>% filter(p.adjust(`pval_trans-PRS`,method="fdr")<0.05) 
outlierscis= cis_wide %>% filter(p.adjust(`pval_cis-PRS`,method="fdr")<0.05)

#### create ranges for axes
range_x = range(c(abs(cis_wide$`cis-z_diffresprot`), abs(trans_wide$`trans-z_diffresprot`)))
range_y = range(c(abs(cis_wide$`zscore_cis-PRS`), abs(trans_wide$`zscore_trans-PRS`)))

###################### SUPPLEMENTARY FIGURE CIS TRANS DIFF-PGS 
diff= ggplot(cis_wide, aes(x = abs(`cis-z_diffresprot`), y = abs(`zscore_cis-PRS`), fill=sig_zdiff)) +
  geom_point(aes(shape=sig_fdr, fill=sig_zdiff), size = 2) + 
  scale_fill_manual(values = c("no" = "lightgrey", "yes" = "blue4"), name = "") + 
  scale_shape_manual(values = c("yes" = 21, "no" = 24), name="FDR < 0.05 in at least\none model") +
  geom_text_repel(data = outlierscis, aes(label = paste0(protein," in\n", Endpoint)),
                  size = 3, vjust = -1, hjust = 1) +
  labs(x = "Difference between cis-adjusted and unadjusted protein ( |z-score| )",
       y = "Association of cis-PGS with incident disease ( |z-score| )") +
  guides(
    fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal() +
  xlim(range_x) + ylim(range_y) +
  theme(
    axis.title.x = element_text(face = "bold",size = 10),
    axis.title.y = element_text(face = "bold",size = 10),
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = 14)
  ) +

  ggplot(trans_wide, aes(x = abs(`trans-z_diffresprot`), y = abs(`zscore_trans-PRS`), fill=sig_zdiff)) +
  geom_point(aes(shape=sig_fdr, fill=sig_zdiff), size = 2) + 
  scale_fill_manual(values = c("no" = "lightgrey", "yes" = "blue4"),
                    name = "Significant difference between\n adjusted and unadjusted protein-disease\n pairs (P < 0.05)") + 
  scale_shape_manual(values = c("yes" = 21, "no" = 24), name="FDR < 0.05 in at least\none model") +
  geom_text_repel(data = outlierstrans, aes(label = paste0(protein," in\n", Endpoint)),
                  size = 3, vjust = -0.5, hjust = 1) +
  labs(x = "Difference between trans-adjusted and unadjusted protein ( |z-score| )",
       y = "Association of trans-PGS with incident disease ( |z-score| )") +
  guides(
    fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal() +
  xlim(range_x) + ylim(range_y) +
  theme(
    axis.title.x = element_text(face = "bold",size = 10),
    axis.title.y = element_text(face = "bold",size = 10),
    plot.tag = element_text(face = "bold", size = 14)
  ) + plot_annotation(tag_levels = "A")

pdf("plots/Supplementary_DIFFPGScistrans.pdf",width = 15, height = 7)
diff
dev.off()

png("plots/Supplementary_DIFFPGScistrans.png",width = 15,height = 7, units = "in",res = 300)
diff
dev.off()

#-------------------------------------------------------------------------------
cis_sig= cis_wide %>% filter(sig_fdr=="yes" & sig_zdiff=="yes")
trans_sig= trans_wide %>% filter(sig_fdr=="yes" & sig_zdiff=="yes")
nrow(cis_sig %>% filter(abs(`zscore_cis-residuals`)>abs(zscore_protein)))
nrow(trans_sig %>% filter(abs(`zscore_trans-residuals`)>abs(zscore_protein)))

nrow(cis_sig %>% filter(`pval_cis-PRS`<0.05))
nrow(trans_sig %>% filter(`pval_trans-PRS`<0.05))


