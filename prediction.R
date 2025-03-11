#' Prediction Functions
#' 
#' Functions for making predictions with the B-spline quantile regression model

#' Create prediction functions based on fitted model
#' 
#' @param sol Solution vector from optimization
#' @param basis_num Number of basis functions in each dimension
#' @param data_dim Dimension of data (number of columns)
#' @param knots_list List of knots for each dimension
#' @param boundary_knots_list List of boundary knots for each dimension
#' @param degree Degree of B-splines
#' @param quantile_func_eval Function to evaluate quantile basis
#' @return List of prediction functions
create_prediction_functions <- function(sol, basis_num, data_dim, knots_list, 
                                        boundary_knots_list, degree, quantile_func_eval) {
  require(splines2)
  
  # Reshape solution into tensor
  sol_tensor <- array(sol, dim = c(rep((degree + 1 + length(knots_list[[1]])), data_dim), 
                                   basis_num[length(basis_num)]))
  
  # Main prediction function
  predict_func <- function(x, p) {
    vec_for_ts_list <- list()
    vec_for_ts_list[[data_dim + 1]] <- quantile_func_eval(p)
    
    for (cnt_pred in 1:data_dim) {
      vec_for_ts_list[[cnt_pred]] <- c(bSpline(
        x[cnt_pred],
        Boundary.knots = boundary_knots_list[[cnt_pred]],
        knots = knots_list[[cnt_pred]],
        degree = degree,
        intercept = TRUE
      ))
    }
    
    # Compute tensor product
    array_tsproduct <- ts_product_func(vec_for_ts_list)
    
    return(sum(sol_tensor * array_tsproduct))
  }
  
  # Vectorized version
  predict_func_vec <- Vectorize(predict_func, vectorize.args = c("x", "p"))
  
  # Sub-function to predict contribution of individual quantile functions
  predict_func_sub <- function(x, index_quantfunc) {
    vec_for_ts_list <- list()
    
    for (cnt_pred in 1:data_dim) {
      vec_for_ts_list[[cnt_pred]] <- c(bSpline(
        x[cnt_pred],
        Boundary.knots = boundary_knots_list[[cnt_pred]],
        knots = knots_list[[cnt_pred]],
        degree = degree,
        intercept = TRUE
      ))
    }
    
    # Compute tensor product
    array_tsproduct <- ts_product_func(vec_for_ts_list)
    
    return(sum(asub(sol_tensor, index_quantfunc, data_dim + 1) * array_tsproduct))
  }
  
  predict_func_sub_vec <- Vectorize(predict_func_sub, vectorize.args = c("x", "index_quantfunc"))
  
  return(list(
    predict_func = predict_func,
    predict_func_vec = predict_func_vec,
    predict_func_sub = predict_func_sub,
    predict_func_sub_vec = predict_func_sub_vec,
    sol_tensor = sol_tensor
  ))
}

#' Generate predictions for test data
#' 
#' @param predict_func_vec Vectorized prediction function
#' @param data_oos Out-of-sample data
#' @param coverage_quantile_vec Vector of quantiles to predict
#' @return Matrix of predicted quantiles
generate_predictions <- function(predict_func_vec, data_oos, coverage_quantile_vec) {
  n_test <- nrow(data_oos)
  n_quantiles <- length(coverage_quantile_vec)
  predictions <- matrix(NA, nrow = n_test, ncol = n_quantiles)
  
  for (i in 1:n_test) {
    predictions[i, ] <- predict_func_vec(
      list(data_oos[i, ]), 
      as.list(coverage_quantile_vec)
    )
  }
  
  return(predictions)
}
