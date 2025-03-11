#' Utility Functions
#'
#' This file contains utility functions used throughout the project.

#' Convert a standard matrix to a sparse matrix format for PSG
#' 
#' @param matrix Input matrix to convert
#' @return Sparse matrix in PSG format
sparse_mat_psg <- function(matrix) {
  sparse_mat <- as(matrix, "dgTMatrix")
  sparse_mat@i = as.integer(1 + sparse_mat@i)
  sparse_mat@j = as.integer(1 + sparse_mat@j)
  return(sparse_mat)
}

#' Function for tensor product of a list of vectors
#' 
#' @param vector_list List of vectors to compute tensor product
#' @return Tensor product result
ts_product_func <- function(vector_list) {
  temp_outer <- vector_list[[1]]
  
  for (cnt in 2:length(vector_list)) {
    # Loop over outer products to create tensor
    temp_outer <- outer(temp_outer, vector_list[[cnt]])
  }
  
  return(temp_outer)
}

#' 
#' 
#' #' Load and preprocess data
#' #' 
#' #' @param file_path Path to the data file
#' #' @return List containing processed data and observations
#' load_preprocess_data <- function(file_path) {
#'   if (!require(readxl)) {
#'     install.packages("readxl")
#'     library(readxl)
#'   }
#'   
#'   data_load <- read_excel(file_path, sheet = 1)
#'   data_load <- as.data.frame(data_load)
#'   observation_all <- data_load[, 5]
#'   data_all <- data_load[, 1:4]
#'   
#'   # Standardize
#'   observation_all <- (observation_all - quantile(observation_all, probs = 0.5)) / 
#'     (quantile(observation_all, probs = 0.75) - quantile(observation_all, probs = 0.25))
#'   
#'   data_all <- apply(data_all, 2, function(x) {
#'     (x - quantile(x, prob = 0.5)) / 
#'       (quantile(x, probs = 0.75) - quantile(x, probs = 0.25))
#'   })
#'   
#'   # Random permutation
#'   random_perm <- sample.int(length(observation_all), length(observation_all))
#'   observation_all <- observation_all[random_perm]
#'   data_all <- as.matrix(data_all[random_perm, ])
#'   
#'   return(list(
#'     observation_all = observation_all,
#'     data_all = data_all
#'   ))
#' }
#' 
#' #' Split data into training and testing sets for cross-validation
#' #' 
#' #' @param observation_all All observations
#' #' @param data_all All data
#' #' @param cnt_cv Current fold number
#' #' @param fold_num Total number of folds
#' #' @return List containing training and testing data
#' create_cv_split <- function(observation_all, data_all, cnt_cv, fold_num) {
#'   # Calculate indices for test set
#'   test_indices <- (1 + (cnt_cv - 1) * length(observation_all) / fold_num):
#'     (cnt_cv * length(observation_all) / fold_num)
#'   
#'   # In-sample data
#'   observation <- observation_all[-test_indices]
#'   data <- as.matrix(data_all[-test_indices, , drop = FALSE])
#'   
#'   # Out-of-sample data
#'   observation_oos <- observation_all[test_indices]
#'   data_oos <- as.matrix(data_all[test_indices, , drop = FALSE])
#'   
#'   return(list(
#'     observation = observation,
#'     data = data,
#'     observation_oos = observation_oos,
#'     data_oos = data_oos
#'   ))
#' }
