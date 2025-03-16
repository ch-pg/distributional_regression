#' Visualization of Distributional Regression
#' 
#' Runs the full analysis including data preprocessing, optimization, prediction, and evaluation.

set.seed(314)

# Source R scripts
source("functions/quantile-functions.R")
source("functions/utils-functions.R")
source("functions/optimization-setup.R")
source("functions/prediction.R")
source("functions/evaluation.R")

# Load libraries
library(readxl)
library(splines2)
library(Matrix)
library(matrixcalc)
library(gurobi)
library(stargazer)
library(matrixStats)
library(plotly)
library(mgcv)

# Initialization ------------------------------------------------------------------------

# Define parameters
# coverage_l <- 0.05
# coverage_u <- 0.95
coverage_quantile_vec <- c(0.01, seq(0.05, 0.95, 0.1), 0.99) # confidence levels to compute coverage
conf_vec <- c(0.01,0.05,0.15,0.35,0.45,0.55,0.65,0.85,0.95,0.99) # confidence levels in discretized CRPS optimization
weight_vec <- c(20, 10, rep(1, 6), 10, 20) # weight of pinball loss for each confidence level
weight_vec <- weight_vec / sum(weight_vec)
penalty_coef <- 1 # penalty coefficient for P-spline
degree <- 3 # degree of basis spline function
num_knots <- 6 # number of splines knots of basis spline function
degree_i <- 2 # NA # degree of spline in basis quantile function, NA represents no spline
num_knots_i <- 2 # number of spline knots of basis quantile function
time_limit <- 300 # time limit for optimization

# Load and preprocess data
data_load <-  as.data.frame(read_excel('data/CPI_data.xlsx', sheet = 1, col_names = c('1', '2', '3')))
observation_all <- data_load[, 3]
data_all <- data_load[, c(1,2)]

# Standardization
observation_all <- (observation_all - quantile(observation_all, probs = 0.5)) /
  (quantile(observation_all, probs = 0.75) - quantile(observation_all, probs = 0.25))
data_all <- apply(data_all, 2, function(x) {
  (x - quantile(x, probs = 0.5)) / (quantile(x, probs = 0.75) - quantile(x, probs = 0.25))
})

# Random permutation
random_perm <- sample.int(length(observation_all), length(observation_all))
observation_all <- observation_all[random_perm]
data_all <- as.matrix(data_all[random_perm, ])

# Initialize storage
n_methods <- 1
crps_oos_mat <- matrix(NA, nrow = n_methods, ncol = 1)
crps_oos_mat_med <- matrix(NA, nrow = n_methods, ncol = 1)
coverage_rate <- matrix(NA, nrow = 1, ncol = length(coverage_quantile_vec))
opt_time <- matrix(NA, nrow = n_methods, ncol = 1)
predict_quantile_list <- list()
for (cnt in 1:n_methods) {
  predict_quantile_list[[cnt]] <- list()
  for (cnt1 in 1:1) {
    predict_quantile_list[[cnt]][[cnt1]] <- matrix(
      NA,
      nrow = length(observation_all) / 1,
      ncol = length(coverage_quantile_vec)
    )
  }
}

# Create quantile functions
quantile_funcs <- create_quantile_func_list(degree_i, num_knots_i, 
                                            quantile_funcs = c("qnorm", 
                                                               "qexp_def", 
                                                               "qexp_def_left"))
quantile_func_list <- quantile_funcs$quantile_func_list
quantile_func_eval <- quantile_funcs$quantile_func_eval
nonspline_num <- quantile_funcs$nonspline_num


# Optimization ------------------------------------------------------------------------

# 
fold_size <- floor(length(observation_all) / 2)
idx_oos <- 1:fold_size

observation <- observation_all[-idx_oos]
data <- as.matrix(data_all[-idx_oos, ])
observation_oos <- observation_all[idx_oos]
data_oos <- as.matrix(data_all[idx_oos, ])

# Create B-spline knots
knots <- create_bspline_knots(data, num_knots)
knots_list <- knots$knots_list
boundary_knots_list <- knots$boundary_knots_list

# Create design matrices
design_mats <- create_design_matrices(
  data, observation, conf_vec, knots_list, boundary_knots_list, degree, quantile_func_eval
)
design_mat_list <- design_mats$design_mat_list
ncol_design_mat <- design_mats$ncol_design_mat
basis_num <- design_mats$basis_num

# Create quadratic penalty matrix
quadratic_mat <- create_quadratic_penalty(basis_num, ncol(data))

# Setup optimization components
setup <- setup_optimization_gurobi(
  design_mat_list, observation, conf_vec, weight_vec, quadratic_mat, 
  num_knots, degree, ncol_design_mat, ncol(data)
)

# Assemble Gurobi model
model <- setup_model_gurobi(
  observation, conf_vec, quadratic_mat, penalty_coef, 
  setup$ncol_constr_mat, setup$quadratic_mat_new, setup$c, 
  setup$constr_mat, setup$lb, setup$ub, time_limit
)

# Optimization parameters
params <- list(
  OptimalityTol = 1e-2,
  Method = -1,
  TimeLimit = time_limit
)

# Solve the model
result <- gurobi(model, params)
sol <- result$x[1:ncol_design_mat]
opt_time[1, 1] <- result$runtime



# Evaluation ------------------------------------------------------------------------

# Create prediction functions
pred_funcs <- create_prediction_functions(
  sol, basis_num, ncol(data), knots_list, boundary_knots_list, degree, quantile_func_eval
)
predict_func <- pred_funcs$predict_func
predict_func_vec <- pred_funcs$predict_func_vec

# Generate predictions
predictions <- generate_predictions(predict_func_vec, data_oos, coverage_quantile_vec)
predict_quantile_list[[1]][[1]] <- predictions

# Compute out-of-sample CRPS
crps_result <- compute_crps(observation_oos, data_oos, predict_func_vec)
crps_oos_mat[1, 1] <- crps_result$mean_crps
crps_oos_mat_med[1, 1] <- crps_result$median_crps

# Compute out-of-sample coverage rates
coverage_rate[1, ] <- compute_all_coverage_rates(
  observation_oos, predict_quantile_list[[1]][[1]], coverage_quantile_vec
)

# # Clean up
# rm(design_mat_list, quadratic_mat)
# gc()


# CRPS
print(crps_oos_mat)

# Plot coverage curves
plot_coverage_curves(coverage_rate, coverage_quantile_vec, mean_cv_coverage = FALSE, na.rm = TRUE)

# Generate 3D plots
generate_3d_plots(data, observation, predict_func, x_feature = 1, y_feature = 2, 
                  alpha = 0.3, confidence_vec_plot = c(0.01, 0.99))