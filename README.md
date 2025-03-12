# Nonparametric Distributional Regression via Quantile Regression

This repository contains an implementation of a nonparametric distributional regression method, along with cross-validation procedures to evaluate its performance.

## Structure

```
distributional_regression/
├── README.md                # Project description
├── data/                    # Data directory
│   └── Folds5x2_pp.xlsx     # Dataset (Combined Cycle Power Plant from UCI Data)
├── functions/               # Functions directory
│   ├── quantile-functions.R # Basis quantile functions definitions
│   ├── utils-functions.R    # Utility functions
│   ├── optimization-setup.R # Variable setup for optimization
│   ├── prediction.R         # Prediction functions
│   └── evaluation.R         # Evaluation metrics and plots (CRPS, coverage rate, 3d plot)
└── main.R                   # Main script that sources and runs everything (cross validation)
```
<!---
├── results/                 # Results directory
│   └── ...                  # Generated plots, tables, etc.
--->

<!---
## Dependencies

The following R packages are required:
- readxl
- Matrix
- splines2
- gurobi
- ggplot2
- reshape2
- abind
- dplyr
--->

## Usage

 <!---
1. Clone the repository
2. Ensure you have Gurobi installed and configured
3. Put your data file in the data directory
4. Run the main script:
--->

```r
source("main.R")
```
 <!---
## Description of Modules

- **utils.R**: Contains utility functions for data manipulation and tensor product operations
- **basis_functions.R**: Functions for creating B-spline basis functions
- **quantile_functions.R**: Functions for creating and evaluating quantile basis functions
- **optimization.R**: Functions for setting up and running the Gurobi optimization
- **prediction.R**: Functions for making predictions with the fitted model
- **evaluation.R**: Functions for evaluating model performance
- **main.R**: Main script that orchestrates the workflow

--->

## Model Description

The proposed model is a new approach to estimating
the distribution of a response variable
conditioned on factors. We model the
conditional quantile function as a mixture
(weighted sum) of basis quantile functions,
with weights depending on these factors.
The estimation problem is formulated as
a convex optimization problem solved efficiently by Gurobi. The objective
function is equivalent to the continuous
ranked probability score (CRPS).
This approach can be viewed as conducting
quantile regressions for all confidence levels
simultaneously while inherently avoiding
quantile crossing. We use spline functions
of factors as a primary example for the
weight function.


<!---
## Evaluation Metrics

The model is evaluated using:
- Pinball loss for each quantile
- Coverage (proportion of true values below predicted quantile)
- Interval widths between consecutive quantiles
- Calibration curves
- CRPS (Continuous Ranked Probability Score)

## License

[Insert your license information here]

--->