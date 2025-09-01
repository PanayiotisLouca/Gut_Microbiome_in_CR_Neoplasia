# Gut Microbiome Dynamics in Colorectal Neoplasia: Associations Across Disease Trajectory and Adenomatous, Serrated, and Sessile Serrated Polyp Subtypes

## Authors 

*Panayiotis Louca, Sarah Manning, Eleanor Hackney, Linda Sharp, Mark A. Hull, Sarah Koo, Greg R. Young, Guy Taylor, Yashvee Dunneram, Suparna Mitra, James Hampton, Christina Dobson, Laura Neilson, Caroline Addison, Emad M El-Omar, the COLOCOHORT research team, Christopher J. Stewart, Colin Rees*

---

  Analysis scripts used in our comprehensive study of the human gut bacteriome across neoplasia categories (DOI: To be updated). 

---

  ## 📂 Repository Structure
```
├── 1.Alpha_diversity_ANOVA
│   └── 1.1.Alpha_Diversity_SCRIPT.R
├── 2.Beta_diversity_PERMANOVA
│   └── 2.1.beta_diversity_SCRIPT.R
├── 3.Random_Forest_Recursive_Feature_Elimination
│   ├── 3.1.RF_feature_reduction_species_repeats.R
│   ├── 3.2.RF_feature_reduction_genus_repeats.R
│   └── 3.3.RF_feature_reduction_phylum_repeats.R
├── 4.Taxonomy_MaAsLin2
│   ├── 4.1.species_Maaslin2_SCRIPT.R
│   ├── 4.2.genus_Maaslin2_SCRIPT.R
│   └── 4.3.phyla_Maaslin2_SCRIPT.R
├── 5.Functional_Potential_MaAsLin2
│   ├── 5.1.gene_families_MaAslin2_SCRIPT.R
│   └── 5.2.pathways_MaAslin2_SCRIPT.R
└── 6.High_risk_RF_classification
    ├── 6.1.high_risk_prediction_SCRIPT_with_feature_importance.R
    └── 6.2.high_risk_prediction_SCRIPT_with_feature_importance_sensitivity_analysis.R
```

  ## ℹ️ Repository Information 

  - **Alpha diversity (ANOVA)**
      - `1.1.Alpha_Diversity_SCRIPT.R`: Assesses differences in alpha diversity metrics (Shannon diversity and observed taxonomic richness) between polyp-negative individuals and those with neoplasia, as well as comparisons across the neoplasia disease spectrum and polyp subtypes.

  - **Beta_diversity (PERMANOVA)**
      - `2.1.beta_diversity_SCRIPT.R`: Uses the `vegan` package and `adonis2` to test differences in gut microbiome community composition (Bray-Curtis and Jaccard dissimilarities).

  - **Random_Forest_Recursive Feature Elimination**
      - `3.1.RF_feature_reduction_species_repeats.R`: Performs recursive feature elimination using random forest with 5-fold cross-validation (100 repeats, different seeds), integrating gut microbial *species*-level data with host covariates.
      - `3.2.RF_feature_reduction_genus_repeats.R`: As above, but for *genus*-level data.
      - `3.3.RF_feature_reduction_phylum_repeats.R`: As above, but for *phylum*-level data.

  - **Taxonomy (MaAsLin2)**
      - `4.1.species_Maaslin2_SCRIPT.R`: Assesses associations between microbial *species* relative abundances and neoplasia outcomes.
      - `4.2.genus_Maaslin2_SCRIPT.R`: As above, for *genera*.
      - `4.3.phyla_Maaslin2_SCRIPT.R`: As above, for *phyla*.

  - **Functional Potential (MaAsLin2)**
      - `5.1.gene_families_MaAslin2_SCRIPT.R`: Uses linear mixed-effects models to test associations between neoplasia and gene families (filtered for UNIPROT proteins with enzyme commission numbers).
      - `5.2.pathways_MaAslin2_SCRIPT.R`: Models associations between neoplasia and MetaCyc pathways.

  - **High-risk RF classification**
      - `6.1.high_risk_prediction_SCRIPT_with_feature_importance.R`: Matches high-risk participants (age-, sex-, and BMI-matched) to polyp-free controls and builds RF models with 5-fold cross-validation (microbiome-only, covariates-only, and combined).
      - `6.2.high_risk_prediction_SCRIPT_with_feature_importance_sensitivity_analysis.R`: Sensitivity analysis using a broader control pool.
          
  ## Citation 
  If you use this code, please cite:
  Louca, P. et al. (2025). *Gut Microbiome Dynamics in Colorectal Neoplasia: Associations Across Disease Trajectory and Adenomatous, Serrated, and Sessile Serrated Polyp Subtypes*. [DOI: To be updated]  
