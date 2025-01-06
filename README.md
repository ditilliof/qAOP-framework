# qAOP Framework

This repository provides resources to generate figures and re-run model calibration for the study:  
*Developing Quantitative Adverse Outcome Pathways: An Ordinary Differential Equation-Based Computational Framework*.  
Published in **Computational Toxicology**. DOI: [10.1016/j.comtox.2024.100330](https://doi.org/10.1016/j.comtox.2024.100330).

## Folder Contents

### `feedback_qAOP`
This folder contains files for generating Figures 1 and 5:
- **Figure 1**: Run the code chunks in `feedback_qAOP.Rmd`. The bifurcation plot (Figure 1B) generation entails the call of `grindnew.R` for phase portrait and bifurcation analyses.  
- **Figure 5**: Run the Python script `sAOP_dist_panel.py`. This script generates Figure 5 and saves individual and composite distribution plots in a timestamped folder within the `simpleAOP` directory.

### `model_comparison`
This folder contains files for generating Figure 2:
- Run the R script `model_comp.R`. This script:
  - Loads the three qAOP model functions.
  - Generates artificial data for the selected model.
  - Prompts the user to choose two models for comparison by entering their corresponding capital letters (e.g., `A`, `B`, or `C`).
  - Prompts the user to select the model from which artificial data should be generated.
  - Performs cmdstanr parameter inference for both models using the generated artificial data.

- Calibration files (posterior draws, parameter summaries, LOO diagnostics, trace plots, and density plots) will be saved in a timestamped folder within the `model_comparison` directory.

- Model calibration is performed using files named `qAOPX.stan`, where `X` corresponds to the model selected (`A`, `B`, or `C`).

- To plot the model fits, run the code chunks in `timeplots.Rmd`.

### `model_updating`
This folder contains files for generating Figure 3:
- Run all the code chunks in `model_updating.Rmd` to:
  - Define qAOP models.
  - Generate artificial data.
  - Run cmdstanr model fitting.

### `response-ODE`
This folder contains files for generating Figure 4:
- Run all the code chunks in `response-ODE.Rmd`.
