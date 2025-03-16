#' Define Quantile Functions
#' 
#' This file contains definitions of various quantile functions
#' used in the B-spline quantile regression model.

#' Exponential quantile function with threshold
#' 
#' @param prob Vector of probabilities
#' @param rate Rate parameter for exponential distribution
#' @return Vector of quantiles
qexp_def <- function(prob, rate = 1) {
  exp_threshold <- 0
  vec <- rep(NA, length(prob))
  vec[which(prob <= exp_threshold)] <- 0
  vec[which(prob > exp_threshold)] <- stats::qexp(
    1 - (1 - prob[which(prob > exp_threshold)]) / (1 - exp_threshold),
    rate)
  return(vec)
}

#' Left-tailed exponential quantile function
#' 
#' @param prob Vector of probabilities
#' @param rate_left Rate parameter for left tail
#' @return Vector of quantiles
qexp_def_left <- function(prob, rate_left = 1) {
  return(-qexp_def(1 - prob, rate_left))
}

#' Lomax quantile function with threshold
#' 
#' @param prob Vector of probabilities
#' @param scale Scale parameter
#' @param shape Shape parameter
#' @return Vector of quantiles
qlomax_def <- function(prob, scale = 1, shape = 3) {
  lomax_threshold <- 0.75
  vec <- rep(NA, length(prob))
  vec[which(prob <= lomax_threshold)] <- 0
  vec[which(prob > lomax_threshold)] <- Renext::qlomax(
    1 - (1 - prob[which(prob > lomax_threshold)]) / (1 - lomax_threshold),
    scale = scale, shape = shape)
  return(vec)
}

#' Left-tailed Lomax quantile function
#' 
#' @param prob Vector of probabilities
#' @param rate_left Rate parameter for left tail
#' @return Vector of quantiles
qlomax_def_left <- function(prob, rate_left = 1) {
  return(-qlomax_def(1 - prob, rate_left))
}




#' Create list of quantile functions
#' 
#' @param degree Degree of I-spline (NA if no splines desired)
#' @param num_knots Number of interior knots for I-spline (NA if no splines desired)
#' @param quantile_funcs List of defined quantile functions of common distributions
#' @return List containing quantile_func_list (list of functions), quantile_func_eval (evaluation function), 
#'         and nonspline_num (number of non-spline functions)
create_quantile_func_list <- function(degree, num_knots, 
                                      quantile_funcs = list(qnorm)) {
  require(splines2)
  
  # Define basis quantile functions using constant function and known distributions
  quantile_func_list <- c(list(function(p) return(1)), lapply(quantile_funcs, match.fun))  
  
  # Verify all elements are functions
  if (!all(sapply(quantile_func_list, is.function))) {
    stop("All provided arguments must correspond to functions")
  }
  
  nonspline_num <- length(quantile_func_list)
  
  # If degree is not NA, prepare for I-spline basis functions
  if (!is.na(degree)) {
    # Placeholder functions for I-spline bases (not actually used, as evaluation is dynamic)
    for (cnt in 1:(num_knots + degree + 1)) {
      quantile_func_list[[nonspline_num + cnt]] <- function(p) return(0)
    }
  }
  
  # Evaluation function
  if (is.na(degree) || is.na(num_knots)) {  # No I-splines
    quantile_func_eval <- function(p) {
      sapply(quantile_func_list, function(func) func(p))[1:nonspline_num]
    }
  } else {  # With I-splines
    quantile_func_eval <- function(p) {
      c(
        sapply(quantile_func_list, function(func) func(p))[1:nonspline_num],
        iSpline(p, Boundary.knots = c(0, 1),
                knots = seq(0, 1, length.out = num_knots + 2)[-c(1, num_knots + 2)],
                degree = degree, intercept = TRUE)
      )
    }
  }
  
  return(list(
    quantile_func_list = quantile_func_list,
    quantile_func_eval = quantile_func_eval,
    nonspline_num = nonspline_num
  ))
}

