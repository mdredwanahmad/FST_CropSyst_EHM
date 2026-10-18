# CropSyst-EHM: Fruit Surface Temperature Modeling and Crop Water Demand Analysis

## Overview

This repository contains the analysis scripts, available datasets, and model files supporting the manuscript:

**“CropSyst-EHM: A Novel Modeling Framework to Quantify Total Crop Water Demand, Including Evapotranspiration and Heat-Stress Management.”**

The materials support fruit surface temperature (FST) model evaluation, meteorological data quality control, irrigation demand estimation, and sensitivity analysis. The repository also includes the CropSyst-EHM executable and crop parameter files used in the simulations.

## Repository contents

```text
.
├── README.md
├── CropSyst_files/
│   ├── Crop/
│   │   ├── apple.crp
│   │   └── blueberry.crp
│   ├── CropSyst_executable/
│   │   └── CropSyst_EHM.exe
│   └── README_source_code.txt
├── data/
│   ├── irrigation_demand/
│   └── [available CSV datasets]
├── DELSA_sensitivity/
│   ├── 01_apple_delsa.ipynb
│   ├── 01_LHS_forFST.ipynb
│   ├── 02_blueberry_delsa.ipynb
│   ├── delsa_test_apple_all.R
│   └── delsa_test_blueberry.R
├── fst_energy_balance_model.R
├── berry_fst_evaluation.R
├── apple_fst_evaluation.R
├── FST_method_comp_AWN.R
├── irrigation_demand_estimation.R
└── qa_qc_sunrise_quincy_corvalis.R
```

### FST modeling and evaluation

| Script | Purpose |
| --- | --- |
| `fst_energy_balance_model.R` | Estimate FST using the fruit energy balance model. |
| `apple_fst_evaluation.R` | Evaluate apple FST simulations, calculate performance statistics, and generate associated figures. |
| `berry_fst_evaluation.R` | Evaluate blueberry and raspberry FST simulations. |
| `FST_method_comp_AWN.R` | Compare FST estimation methods using AgWeatherNet (AWN) meteorological data. |
| `qa_qc_sunrise_quincy_corvalis.R` | Perform quality assurance and quality control (QA/QC) of observational FST data, including screening extreme thermocouple readings. |
| `irrigation_demand_estimation.R` | Analyze irrigation demand, crop phenology, and associated model outputs. |

### Data (`data/`)

The `data/` directory contains available FST observations, meteorological inputs, and model outputs used in the analyses. The `data/irrigation_demand/` subdirectory contains data for the irrigation demand analysis.

**Data availability and limitations**

- Available apple FST observations and associated datasets are included.
- The blueberry observational dataset is **not redistributed** here. Obtain the relevant data from Yang et al. (2024) before running analyses that depend on it.
- Individual AWN station datasets are **not included** because of their size. Contact the corresponding author of Savalkar et al. (2024) for access to the datasets used in that study. However, the processed dataset shared in directory and Appendix C02 figure can be created.

Some analyses therefore require external data in addition to the files provided here. Please consult the manuscript for the complete data sources and references.

## Running the analyses

1. Download or clone this repository, keeping the folder structure unchanged.
2. Install the required R and Python packages identified by the `library(...)` statements and imports in the relevant scripts and notebooks.
3. Open `fst_energy_balance_model.R` script in Rstudio and run this script first to run other scripts.
3. Update your local working directories and file paths in the `apple_fst_evaluation.R` script to point to your copy of the repository. 
4. Run the relevant scripts or notebooks in the order described above.
5. Confirm that all required datasets are present, including externally obtained datasets where applicable.

This section will generate Figure 5 , Figure 8 , and all appendix plots.

### DELSA sensitivity analysis (`DELSA_sensitivity/`)

This folder contains the R scripts and Jupyter notebooks used for parameter-space sampling, model bias calculations, and Distributed Evaluation of Local Sensitivity Analysis (DELSA).

| File | Purpose |
| --- | --- |
| `01_LHS_forFST.ipynb` | Generate FST parameter combinations using Latin Hypercube Sampling (LHS). |
| `delsa_test_apple_all.R` | Estimate apple FST model bias across sampled parameter combinations. |
| `01_apple_delsa.ipynb` | Perform DELSA sensitivity analysis for apple and produce Figure 6. |
| `delsa_test_blueberry.R` | Estimate blueberry FST model bias across sampled parameter combinations. |
| `02_blueberry_delsa.ipynb` | Perform DELSA sensitivity analysis for blueberry and produce Figure 7. |

**Suggested workflow for Figure 6 (apple)**

1. Run the relevant final section of `apple_fst_evaluation.R` to prepare FST data for bias estimation.
2. Run `01_LHS_forFST.ipynb` to generate the sampled parameter sets.
3. Run `delsa_test_apple_all.R` to estimate bias for the sampled parameter sets.
4. Run `01_apple_delsa.ipynb` for the DELSA analysis and Figure 6.

**Suggested workflow for Figure 7 (blueberry)**

1. Obtain the required blueberry dataset and run the corresponding portion of `berry_fst_evaluation.R` to prepare FST data for bias estimation.
2. Run `01_LHS_forFST.ipynb` to generate the sampled parameter sets.
3. Run `delsa_test_blueberry.R` to estimate bias for the sampled parameter sets.
4. Run `02_blueberry_delsa.ipynb` for the DELSA analysis and Figure 7.

> **Note:** Verify the required inputs, file names, and output paths in each script or notebook before running it. Some workflows depend on outputs generated by earlier steps.

### CropSyst model files (`CropSyst_files/`)

- `Crop/apple.crp`: crop parameter file for apple simulations.
- `Crop/blueberry.crp`: crop parameter file for blueberry simulations.
- `CropSyst_executable/CropSyst_EHM.exe`: Windows executable used for the CropSyst-EHM simulations.
- `README_source_code.txt`: additional information on the source code and where it can be accessed.

The corresponding CropSyst source code is available from an SVN repository. See `CropSyst_files/README_source_code.txt` for the repository location and add the exact SVN revision used to compile this executable, if available.


## Citation

If you use materials from this repository, please cite the associated manuscript:

*CropSyst-EHM: A Novel Modeling Framework to Quantify Total Crop Water Demand, Including Evapotranspiration and Heat-Stress Management.

## Contact

For questions about the model, data, or analysis scripts, please contact the corresponding author of the manuscript.
