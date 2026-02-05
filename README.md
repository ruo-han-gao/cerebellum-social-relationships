# Cerebellum & social relationships

This project examined how quality of social relationships affects cerebellar morphology and cerebellar-cerebral covariance in development and aging.

## Steps for preprocessing and analyses

### 1. Cerebellar parcellation

In `script/1_parcellation/`:

- `1_fsl_func_warp_fusion.sh` warps MNI-aligned atlas (Nettekoven et al., 2024) to subject space (via FSL)
- `2_func_native_vols.py` counts volumes per parcel, including outlier detection.

Expected outputs: Dataframe and images of subject-space and functional parcellations.

Expected runtime: ~5 minutes per subject.

### 2. Data preparation

In `script/2_prep/`:

- `1_qc_bad_euler.R` performs QC based on the Euler index of cerebral cortex surface reconstructions.
- `2_hcp_data_prep.R` combines all available quality brain imaging data with social relationship measures into a single analysis-ready data frame.

### 3. Factor analysis

In `script/3_factor_analysis/`:

- `1_hcp_factor_analysis.Rmd` checks the factoribility of the data, conducts exploritory factor analysis, and compares the factor scores between groups.

### 4. General additive modelling

In `scripts/4_gams/`:

- `1_gam_functions.R` defines helper functions for generalized additive model (GAM) fitting and result extraction.
- `2_gam_for_SR.R` sources `1_gam_functions.R` and fits hierarchical GAMs to evaluate the contribution of each independent variable, assess concurvity, and extract partial R and derivatives.

### 5. Analysis of cerebello-cerebral covariance

In `scripts/5_covariance/`:

- `1_covar_analysis.R` checks assumptions and runs linear regression.
- `2_covar_visualization.R` viisualizes covariance effects using heatmaps and line plots.
