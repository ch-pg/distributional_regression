#' Main Script for Distributional Regression
#' 
#' Runs the full analysis including data preprocessing, cross-validation, optimization, prediction, and evaluation.

set.seed(314)

# Source R scripts
source("/Users/cheng/Desktop/FMQ/code/Github/src/quantile-functions.R")
source("/Users/cheng/Desktop/FMQ/code/Github/src/utils-functions.R")
source("/Users/cheng/Desktop/FMQ/code/Github/src/optimization-setup.R")
source("/Users/cheng/Desktop/FMQ/code/Github/src/prediction.R")
source("/Users/cheng/Desktop/FMQ/code/Github/src/evaluation.R")

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

# Define parameters
fold_num <- 3 # number of folds in cross validation
# coverage_l <- 0.05
# coverage_u <- 0.95
coverage_quantile_vec <- c(0.01, seq(0.05, 0.95, 0.1), 0.99) # confidence levels to compute coverage
conf_vec <- c(0.01, 0.33, 0.66, 0.99) # confidence levels in discretized CRPS optimization
weight_vec <- 0.1 * c(5, 2, 2, 5) # weight of pinball loss for each confidence level
penalty_coef <- 0.01 # 0.1 # penalty coefficient for P-spline
degree <- 2 # degree of basis spline function
num_knots <- 1 # number of splines knots of basis spline function
degree_i <- NA # degree of spline in basis quantile function, NA represents no spline
num_knots_i <- 1 # number of spline knots of basis quantile function
time_limit <- 300 # time limit for optimization

# Load and preprocess data
data_load <- read_excel('/Users/cheng/Desktop/FMQ/code/CCPP/Folds5x2_pp.xlsx', sheet = 1)
data_load <- as.data.frame(data_load)
observation_all <- data_load[, 5]
data_all <- data_load[, 1:4]

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
crps_oos_mat <- matrix(NA, nrow = n_methods, ncol = fold_num)
crps_oos_mat_med <- matrix(NA, nrow = n_methods, ncol = fold_num)
coverage_rate <- matrix(NA, nrow = fold_num, ncol = length(coverage_quantile_vec))
opt_time <- matrix(NA, nrow = n_methods, ncol = fold_num)
predict_quantile_list <- list()
for (cnt in 1:n_methods) {
  predict_quantile_list[[cnt]] <- list()
  for (cnt1 in 1:fold_num) {
    predict_quantile_list[[cnt]][[cnt1]] <- matrix(
      NA,
      nrow = length(observation_all) / fold_num,
      ncol = length(coverage_quantile_vec)
    )
  }
}

# Create quantile functions
quantile_funcs <- create_quantile_func_list(degree_i, num_knots_i)
quantile_func_list <- quantile_funcs$quantile_func_list
quantile_func_eval <- quantile_funcs$quantile_func_eval
nonspline_num <- quantile_funcs$nonspline_num

# Cross-validation loop
fold_size <- floor(length(observation_all) / fold_num)

for (cnt_cv in 1:fold_num) {
  cat(sprintf("Fold %d of %d\n", cnt_cv, fold_num))
  
  # Data splitting
  if (cnt_cv < fold_num) {
    # Regular folds
    idx_oos <- (1 + (cnt_cv - 1) * fold_size):(cnt_cv * fold_size)
  } else {
    # Last fold: take remaining observations
    idx_oos <- (1 + (cnt_cv - 1) * fold_size):length(observation_all)
  }
  
  observation <- observation_all[-idx_oos]
  data <- as.matrix(data_all[-idx_oos, ])
  observation_oos <- observation_all[idx_oos]
  data_oos <- as.matrix(data_all[idx_oos, ])
  
  # Create B-spline knots
  knots <- create_bspline_knots(data_all, num_knots)
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
  opt_time[1, cnt_cv] <- result$runtime
  
  # Create prediction functions
  pred_funcs <- create_prediction_functions(
    sol, basis_num, ncol(data), knots_list, boundary_knots_list, degree, quantile_func_eval
  )
  predict_func <- pred_funcs$predict_func
  predict_func_vec <- pred_funcs$predict_func_vec
  
  # Generate predictions
  predictions <- generate_predictions(predict_func_vec, data_oos, coverage_quantile_vec)
  predict_quantile_list[[1]][[cnt_cv]] <- predictions
  
  # Compute CRPS
  crps_result <- compute_crps(observation_oos, data_oos, predict_func_vec)
  crps_oos_mat[1, cnt_cv] <- crps_result$mean_crps
  crps_oos_mat_med[1, cnt_cv] <- crps_result$median_crps
  
  # Compute coverage rates
  coverage_rate[cnt_cv, ] <- compute_all_coverage_rates(
    observation_oos, predict_quantile_list[[1]][[cnt_cv]], coverage_quantile_vec
  )
  
  # # Clean up
  # rm(design_mat_list, quadratic_mat)
  # gc()
}

# CRPS
print(rowMeans2(crps_oos_mat))

# Plot coverage curves
plot_coverage_curves(coverage_rate, coverage_quantile_vec, mean_cv_coverage = FALSE, na.rm = TRUE)

# Generate 3D plots (using last fold's model)
generate_3d_plots(data, observation, predict_func, x_feature = 1, y_feature = 4)