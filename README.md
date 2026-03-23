
# README: Spatially-explicit predictions of future fire occurrence in Victoria, Australia

**Scripts and data for analysis of the manuscript:** *Spatially-explicit predictions of future fire occurrence in Victoria, Australia*

**Authors:** Charlie Hart, Andrew Dowdy, Sarah C. McColl-Gausden, Hamish Clarke, Luke Collins, Amelia French, Trent D. Penman, Nevil Amos, Angie Haslem, Cindy E. Hauser, Jim Thomson, Josephine MacHunter, Matt White & William L. Geary

**Corresponding author:** William L. Geary (billy.geary@unimelb.edu.au)

---

## Description of the Data and File Structure
This repository contains the code and key data required to replicate the modelling and figures presented in the manuscript.

**Note on Replication:** While the full workflow is documented for transparency, Scripts 1–4 and 9 rely on external raw data (DEECA spatial layers, AGCD/AWO climate grids, etc.) that are not included in this repository due to size and licensing restrictions.

To replicate the core results:
* **Run Step 5** to fit the Boosted Regression Trees (BRTs).
* **Run Step 6** to prepare cross-validation data.
* **Run Step 8** to consolidate and summarise the cross-validation performance.
* **Run Step 10** to generate the spatial fire probability predictions.
* **Run Step 11** to recreate the final publication figures.

---

## Folder: scripts
This folder contains the sequential workflow for the analysis.

* **`fun_fit_forced_brt.R`** - Helper functions for fitting BRT models.
* **`Step 1 Make analysis mask.R`** - Generates the 75m study region mask and stratified random sampling points.
* **`Step 2_Make burned area rasters.R`** - Processes Victorian Fire History into annual binary burned area layers.
* **`Step 3a_SPEI from AWO grids.R`** - Calculates 12 and 24-month SPEI from climate data.
* **`Step 3b_Make covariates.R`** - Processes fire weather, ignitions, topography, and Time Since Fire (TSF) predictors.
* **`Step 3c_Make covariate stacks.R`** - Harmonises all spatial layers into unified annual stacks.
* **`Step 3d_Covariate VIc.R`** - Conducts Variance Inflation Factor (VIF) analysis.
* **`Step 4a_Extract data.R`** - Extracts covariate values to sampling points for model training.
* **`Step 5_Fit the model_lr_0.005.R`** - Fits primary BRT models (learning rate 0.005).
* **`Step 6_Prep Cross validation Data.R`** - Prepares data folds for model validation.
* **`Steps 7a, 7b, 7c (Cross Validation)`** - Batch submission and worker scripts to run model cross-validation in parallel on Spartan.
* **`Step 8_Summarise CVs.R`** - Consolidates performance metrics from the cross-validation runs.
* **`Steps 9a, 9b, 9c (Future Climate Prep)`** - Downscales climate projections, generates future SPEI maps, and compiles future prediction stacks.
* **`Step 10_Current and Future Fire Predictions.R`** - Generates spatial probability rasters for baseline and future scenarios.
* **`Step 11_Plots for publication.R`** - Produces final binned maps, delta maps, and regional boxplots.

**Note on additional_scripts:** This subfolder contains ad-hoc workflows not required to recreate the primary outputs. These include scripts for downloading raw AWO climate data, rolling cross-validation routines, and preliminary SPEI raster generation. They are included for additional context and methodological completeness.

---

## Folder: data
This folder contains essential input files for the analysis.
* **Modelling Data:** Extracted covariate values used to train the BRT models.
* **Districts:** Spatial boundaries for Victorian fire management districts.
* **Random Points:** The stratified random points used for data extraction and model training.

---

## Folder: outputs
This folder contains final datasets and model objects created at key stages.

* **Prediction_Stacks:** Multi-layer annual covariate stacks used to generate predictions.
* **Predictions:** Baseline and future fire probability GeoTIFFs (`.tif`).
* **Final_Models:** * `brt_base_lr0.005_ffdi_95.RDS` - Base model object.
    * `brt_full_rac_lr0.005_ffdi_95.rds` - RAC model object.
* **Regional_Summaries:** * `regional_fire_risk_summary_RAC.csv` - Mean fire probabilities by district.

---

## Code and Software
All analyses were performed in **R v4.3.0** or later. Key libraries: `terra`, `sf`, `ggplot2`, `tidyterra`, `dismo`, and `gbm`.

