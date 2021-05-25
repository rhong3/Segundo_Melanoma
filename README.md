# Melanoma Multi-omics Data Analyses

## Scripts
### Data_prep
 - Cleaning.R: preprocessing for proteomics, phosphoproteomics, transcriptomics, clinical, and histological data
 - ICA_subtype_prep.R: prepare data frame to add proteomics subtypes for ICA
 - post_ICA.R: ICA output formatting  
 - OLA_prep.R: prepare data frame for outlier analyses
 - IHC_model.R: prepare IHC validation data frames; build IHC validation regression models; plot IHC validation results

### Figures
 - Aggregate.R: summarize results from ICA, outlier analyses, and Cox regression analyses
 - density.R: make density plots of proteomics and phosphoproteomics data for each sample
 - graph.R: make graphs to represent the relationship between ICs and pathways
 - graphlite.R: make graphs to represent the relationship between ICs and pathways (filtered)
 - GSEA.R: IC pre-ranked GSEA; make heatmap of the top-ranked proteins
 - ica_pca_comparison.R: ICA and PCA analyses comparison
 - Kaplan-Meier6.R: draw KM curves based on different AUROC cutoffs
 - OLA_HM.R: make heatmap showing the expression level of all the outliers found by outlier analyses across all samples
 - QC_figure.R: quality control for omics data (Figure S1)
 - RNA-PROT-cor.R: protein and RNA expression correlations
 - top-to-pathway.R: significant pathways found by ICA that involve the 10 selected proteins

### ICA
 - clinical_association_template.R: find association between ICs and clinical variables
 - ICA_Clusters_Functions.R: ICA clustering helper functions
 - ICA_subtype_HM.R: Heatmaps of top proteins' expression ranked by ICs that correlate with the proteomics subtypes
 - ICA_template.R: independent component analysis

### outlier_analysis
 - outlier analysis (Blacksheep); please refer to the document files in the folder
 - https://www.biorxiv.org/content/10.1101/825067v2.full.pdf

