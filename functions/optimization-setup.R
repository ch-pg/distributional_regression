#' Optimization Setup Functions for Gurobi
#' 
#' Functions for preparing inputs and setting up the Gurobi optimization problem for distributional regression.

#' Create a permutation matrix function
#' 
#' @param vec_old Original vector
#' @param vec_new New vector with permuted order
#' @return Permutation mapping (indices mapping vec_old to vec_new)
permu_mat_func <- function(vec_old, vec_new) {
  return(sapply(vec_old, function(x) which(vec_new == x)))
}

#' Create knots for B-spline basis
#' 
#' @param data Input data matrix
#' @param num_knots Number of interior knots
#' @return List containing knots_list (interior knots) and boundary_knots_list (boundary knots)
create_bspline_knots <- function(data, num_knots) {
  knots_list <- list()
  boundary_knots_list <- list()
  
  for (cnt in 1:ncol(data)) {
    knots_list[[cnt]] <- seq(
      min(data[, cnt]), 
      max(data[, cnt]), 
      length.out = num_knots + 2
    )[-c(1, num_knots + 2)]
    boundary_knots_list[[cnt]] <- range(data[, cnt])
  }
  
  return(list(
    knots_list = knots_list,
    boundary_knots_list = boundary_knots_list
  ))
}

#' Create design matrices for quantile regression
#' 
#' @param data Input data matrix
#' @param observation Vector of observed values
#' @param conf_vec Vector of confidence levels
#' @param knots_list List of interior knots for each variable
#' @param boundary_knots_list List of boundary knots for each variable
#' @param degree Degree of B-spline
#' @param quantile_func_eval Function to evaluate quantile basis
#' @return List containing design_mat_list (design matrices), ncol_design_mat (number of columns), and basis_num (basis sizes)
create_design_matrices <- function(data, observation, conf_vec, 
                                   knots_list, boundary_knots_list, 
                                   degree, quantile_func_eval) {
  require(splines2)
  
  basis_num <- c(
    unlist(lapply(knots_list, length)) + degree + 1, 
    length(eval(parse(text = "quantile_func_list")))
  )
  ncol_design_mat <- prod(basis_num)
  
  design_mat_pre <- matrix(NA, nrow = length(observation), ncol = ncol_design_mat)
  design_mat_list <- list()
  
  for (cnt in 1:length(conf_vec)) {
    for (cnt0 in 1:nrow(data)) {
      vec_for_ts_list <- list()
      vec_for_ts_list[[ncol(data) + 1]] <- quantile_func_eval(conf_vec[cnt])
      for (cnt1 in 1:ncol(data)) {
        vec_for_ts_list[[cnt1]] <- c(bSpline(
          data[cnt0, cnt1],
          Boundary.knots = boundary_knots_list[[cnt1]],
          knots = knots_list[[cnt1]],
          degree = degree,
          intercept = TRUE
        ))
      }
      array_tsproduct <- ts_product_func(vec_for_ts_list)
      design_mat_pre[cnt0, ] <- c(array_tsproduct)
    }
    design_mat_list[[cnt]] <- cbind(design_mat_pre, observation)
  }
  
  return(list(
    design_mat_list = design_mat_list,
    ncol_design_mat = ncol_design_mat,
    basis_num = basis_num
  ))
}

#' Create quadratic penalty matrix for P-splines
#' 
#' @param basis_num Vector with number of basis functions in each dimension
#' @param ncol_data Number of columns in data
#' @return Quadratic penalty matrix
create_quadratic_penalty <- function(basis_num, ncol_data) {
  require(Matrix)
  require(matrixcalc)
  require(mgcv)  # For sdiag
  
  first_axis <- 1
  quadratic_sub_sub_mat <- matrix(0, basis_num[-length(basis_num)][first_axis], 
                                  basis_num[-length(basis_num)][first_axis])
  
  sdiag(quadratic_sub_sub_mat, 0) <- c(1, 5, rep(6, ncol(quadratic_sub_sub_mat) - 4), 5, 1)
  sdiag(quadratic_sub_sub_mat, -1) <- c(-4, rep(-8, ncol(quadratic_sub_sub_mat) - 3), -4)
  sdiag(quadratic_sub_sub_mat, -2) <- rep(2, ncol(quadratic_sub_sub_mat) - 2)
  
  quadratic_sub_mat <- bdiag(
    rep(list(quadratic_sub_sub_mat), prod(basis_num[-length(basis_num)][-first_axis]))
  )
  
  for (cnt in 2:ncol_data) {
    perm_temp <- c(1:ncol_data)
    perm_temp[c(1, cnt)] <- perm_temp[c(cnt, 1)]
    array_to_perm <- array(
      1:prod(basis_num[-length(basis_num)]), 
      dim = c(basis_num[-length(basis_num)])
    )
    array_permed <- aperm(array_to_perm, perm_temp)
    perm_vec <- permu_mat_func(
      c(1:prod(basis_num[-length(basis_num)])), 
      c(array_permed)
    )
    quadratic_sub_mat <- quadratic_sub_mat + quadratic_sub_mat[perm_vec, perm_vec]
  }
  
  quadratic_sub_mat <- (quadratic_sub_mat + t(quadratic_sub_mat)) / 2
  quadratic_mat <- bdiag(rep(list(quadratic_sub_mat), basis_num[length(basis_num)]))
  
  return(quadratic_mat)
}

#' Set up Gurobi optimization problem components
#' 
#' @param design_mat_list List of design matrices
#' @param observation Vector of observed values
#' @param conf_vec Vector of confidence levels
#' @param weight_vec Vector of weights for different confidence levels
#' @param quadratic_mat Quadratic penalty matrix
#' @param num_knots Number of interior knots
#' @param degree Degree of B-splines
#' @param ncol_design_mat Number of columns in design matrix
#' @param ncol_data Number of columns in data
#' @return List containing optimization components (ncol_constr_mat, quadratic_mat_new, c, constr_mat, lb, ub)
setup_optimization_gurobi <- function(design_mat_list, observation, conf_vec, weight_vec,
                                      quadratic_mat, num_knots, degree, 
                                      ncol_design_mat, ncol_data) {
  ncol_constr_mat <- ncol_design_mat + 2 * length(conf_vec) * length(observation)
  
  zero_matrix <- Matrix(0, 
                        nrow = 2 * length(conf_vec) * length(observation),
                        ncol = 2 * length(conf_vec) * length(observation), 
                        sparse = TRUE)
  quadratic_mat_new <- bdiag(quadratic_mat, zero_matrix)
  
  c <- rep(0, ncol_constr_mat)
  for (k in 1:length(conf_vec)) {
    idx_u <- ncol_design_mat + (k - 1) * 2 * length(observation) + 1
    idx_v <- idx_u + length(observation)
    c[idx_u:(idx_u + length(observation) - 1)] <- conf_vec[k] * weight_vec[k]
    c[idx_v:(idx_v + length(observation) - 1)] <- (1 - conf_vec[k]) * weight_vec[k]
  }
  
  constr_mat_1 <- do.call(rbind, design_mat_list)
  constr_mat_1 <- constr_mat_1[, -ncol(constr_mat_1)]
  identity_mat <- sparseMatrix(
    i = 1:length(observation), 
    j = 1:length(observation), 
    x = 1
  )
  identity_mat_cb <- cbind(identity_mat, -identity_mat)
  constr_mat_2 <- bdiag(replicate(length(conf_vec), identity_mat_cb, simplify = FALSE))
  constr_mat <- cbind(constr_mat_1, constr_mat_2)
  
  lb <- c(
    rep(-Inf, (num_knots + degree + 1)^ncol_data),
    rep(0, ncol_design_mat - (num_knots + degree + 1)^ncol_data),
    rep(0, 2 * length(observation) * length(conf_vec))
  )
  ub <- rep(Inf, ncol_constr_mat)
  
  return(list(
    ncol_constr_mat = ncol_constr_mat,
    quadratic_mat_new = quadratic_mat_new,
    c = c,
    constr_mat = constr_mat,
    lb = lb,
    ub = ub
  ))
}

#' Assemble Gurobi model for quantile regression
#' 
#' @param observation Vector of observed values
#' @param conf_vec Vector of confidence levels
#' @param quadratic_mat Quadratic penalty matrix
#' @param penalty_coef Coefficient for quadratic penalty
#' @param ncol_constr_mat Number of columns in constraint matrix
#' @param quadratic_mat_new Augmented quadratic matrix with slacks
#' @param c Linear objective coefficients
#' @param constr_mat Constraint matrix
#' @param lb Lower bounds for variables
#' @param ub Upper bounds for variables
#' @param time_limit Time limit for optimization in seconds (default: 300)
#' @return Gurobi model list
setup_model_gurobi <- function(observation, conf_vec, quadratic_mat, penalty_coef, 
                               ncol_constr_mat, quadratic_mat_new, c, constr_mat, 
                               lb, ub, time_limit = 300) {
  model <- list()
  model$Q <- quadratic_mat_new * penalty_coef
  model$obj <- c
  model$A <- constr_mat
  model$rhs <- rep(observation, length(conf_vec))
  model$sense <- rep('=', nrow(constr_mat))
  model$lb <- lb
  model$ub <- ub
  model$vtype <- 'C'
  model$modelsense <- 'min'
  
  rm(quadratic_mat_new, constr_mat)
  gc()
  
  return(model)
}