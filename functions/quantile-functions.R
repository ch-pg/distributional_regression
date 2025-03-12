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

#' Create list of basis quantile functions
#' 
#' @param degree_i Degree of spline basis (NA for no spline basis)
#' @param num_knots_i Number of interior knots
#' @return List containing quantile function evaluation function
create_quantile_func_list <- function(degree_i = NA, num_knots_i = 1) {
  # Basic quantile functions (variable name not changed), only right tail
  quantile_func_list <- list(function(p) return(1), qnorm)
  quantile_func_list_defined_length <- length(quantile_func_list)
  nonspline_num <- length(quantile_func_list)
  
  # Add placeholder functions if using splines
  if (!is.na(degree_i)) {
    for (cnt in 1:(num_knots_i + degree_i + 1)) {
      quantile_func_list[[cnt + quantile_func_list_defined_length]] <- function(p) {
        return(0)
      } # place holder
    }
  }
  
  # Helper function to evaluate with arguments
  eval.with.args <- function(FUN, ...) FUN(...)
  
  # Create evaluation function
  if (is.na(degree_i + num_knots_i)) { # no I-spline
    quantile_func_eval <- function(p) {
      return(unlist(lapply(quantile_func_list, eval.with.args, p))[1:nonspline_num])
    }
  } else {
    quantile_func_eval <- function(p) {
      return(c(
        unlist(lapply(quantile_func_list, eval.with.args, p))[1:nonspline_num],
        iSpline(p, Boundary.knots = c(0, 1),
                knots = seq(0, 1, length.out = num_knots_i + 2)[-c(1, num_knots_i + 2)],
                degree = degree_i, intercept = TRUE)
      ))
    }
  }
  
  # Return a list with both components
  return(list(
    quantile_func_list = quantile_func_list,
    quantile_func_eval = quantile_func_eval,
    nonspline_num = nonspline_num
  ))
}
