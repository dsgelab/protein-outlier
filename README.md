# Removing genetic effects on plasma proteins enhances their utility as predictive disease biomarkers

This repository contains code and analyses for the manuscript *"Removing genetic effects on plasma proteins enhances their utility as predictive disease biomarkers"*.

## Data
Data used in this study are available under request via [UK Biobank](https://www.ukbiobank.ac.uk/) and [FinnGen](https://www.finngen.fi/).

## Repository Structure

### `genetic_adjustment/`

Scripts for generating genetically adjusted protein levels.

- `snakefile_pgs_olink` and `snakefile_pgs_somalogic` - Code for generating pPGS from [OmicsPred](https://www.omicspred.org/) Olink and Somalogic derived genetic scores.

- `proteins_olinksoma_filtering.R` (that uses `R2_pPGSprot.R`) - Code to filter proteins based on R<sup>2</sup> with pPGS and their correlations. 

- `snakefile_residualizing` (that uses `create_inputfiles_for_residuals.R` and `ProteinPgsResidual.cpp`) - Code for creating "genetically adjusted proteins", here referred to as "residuals".

- `snakefile_pgs_cistrans` - Code for generating cis- and trans- pPGS and genetically adjusted proteins based on loci defined here `proteins_loci.R`.

### `disease_association/`

Scripts to run association analyses with disease endpoints.

- `incident_logistic_singleprot.R`- Code for incident disease analysis results as in **Figure 2**, **Supplementary Figures 3** and **5** and **Supplementary Table 2**.

- `power-calculation.R` - Code to generate **Supplementary Table 3**.

- `prevalent_logistic_singleprot.R`- Code for prevalent cases analysis results as in **Supplementary Figure 4** and **Supplementary Table 4**.

- `cistrans_logistic_singleprot.R` - Code for cis- and trans- adjusted proteins incident disease analysis results as in **Supplementary Figures 7** and **8** and **Supplementary Table 5**.

### `sensitivity_analysis_filtered_pPGS/`

Scripts to perform sensitivity analysis excluding disease-associated variants from pPGS.

- `snakefile_filtering_pgs`- Code to filter OmicsPred files, compute pPGS and genetically adjusted proteins with filtered SNPs.

- `revision_instrument_selection.R` - Code to produce **Supplementary Table 6**.

- `V2_revision_inc_logreg_afterinstrumentselection.R` - Code for incident disease analysis results with filtered adjusted proteins and pPGS, as in **Supplementary Figure 11** and **Supplementary Table 7**.

### `exposome_analysis/`

Scripts to perform associations analyses with the exposome. 

- `enviromental_exposures_regressions.R`- Code for associations between unadjusted and adjusted proteins with the exposome as in **Figure 4** and **Supplementary Table 8**.

- `revision_cis_environmental_exposures_regressions.R`- Code for associations between cis-adjusted proteins with the exposome as in **Supplementary Figure 12** and **Supplementary Table 8**.

### `multi_protein_models/`

Scripts to perform LASSO-based multi-protein analyses.

- `LASSO_AUC.R`- Code for multi-protein AUC extraction and calibration as in **Figure 5**, **Supplementary Figures 14** and **15**, **Supplementary Table 9**.

- `LASSO_coxmodel.R`- Code for multi-protein C-index extraction and calibration as in **Supplementary Figure 13** and **Supplementary Table 10**.

### `replication_with_finngen_derived_pPGSweights/`

Scripts to perform replication using independently derived pPGS weights from FinnGen including 7 additional proteins.

- `remove_correlation_FGproteins.R` - Code to filter out proteins correlated with the original set of 94 proteins.

- `FGproteins_PGS_LDpred.R`and `LDpredest_per_chr_func.R` - Code to extract weights with LDpred from FinnGen proteins' summary statistics. 

- `snakefile_replication_pgs `- Code to perform pPGS and genetically adjusted proteins in the UKB using FinnGen derived weights.

- `filtering_FG_proteins.R`- Code to filter proteins based on R<sup>2</sup> with pPGS

- `FG_incident_logistic_singleprot.R` - Code for incident disease analysis results as in **Figure 6**, **Supplementary Figure 16** and **Supplementary Table 11**.







