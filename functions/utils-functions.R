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

