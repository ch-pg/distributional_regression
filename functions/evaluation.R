#' Evaluation Functions
#' 
#' Functions for evaluating the B-spline quantile regression model, including CRPS, coverage rates, and visualization.

#' Compute CRPS for out-of-sample predictions
#' 
#' @param observation_oos Out-of-sample observations (vector)
#' @param data_oos Out-of-sample data (matrix)
#' @param predict_func_vec Vectorized prediction function
#' @return List containing mean_crps, median_crps, and crps_vec
compute_crps <- function(observation_oos, data_oos, predict_func_vec) {
  epsilon_integ <- 0.0001
  
  integrand_func <- function(p, x, y) {
    prediction <- predict_func_vec(list(x), list(p))
    if (y - prediction >= 0) {
      return(p * (y - prediction))
    } else {
      return(-(1 - p) * (y - prediction))
    }
  }
  
  crps_vec <- rep(NA, length(observation_oos))
  
  for (cnt in 1:length(observation_oos)) {
    integrand_func_temp_vec <- Vectorize(
      function(p) integrand_func(p, x = data_oos[cnt, ], y = observation_oos[cnt]),
      vectorize.args = "p"
    )
    crps_vec[cnt] <- 2 * integrate(
      integrand_func_temp_vec,
      lower = epsilon_integ,
      upper = 1 - epsilon_integ
    )$value
  }
  
  return(list(
    mean_crps = mean(crps_vec),
    median_crps = median(crps_vec),
    crps_vec = crps_vec
  ))
}

#' Compute coverage rates for prediction intervals
#'
#' @param observation_oos Out-of-sample observations (vector)
#' @param data_oos Out-of-sample data (matrix)
#' @param predict_func_vec Vectorized prediction function
#' @param coverage_l Lower coverage bound (scalar)
#' @param coverage_u Upper coverage bound (scalar)
#' @return Coverage rate (scalar)
compute_coverage_rate <- function(observation_oos, data_oos, predict_func_vec, coverage_l, coverage_u) {
  coverage_binary <- rep(NA, length(observation_oos))

  for (cnt in 1:length(observation_oos)) {
    predict_interval <- predict_func_vec(list(data_oos[cnt, ]), list(coverage_l, coverage_u))
    if (observation_oos[cnt] < predict_interval[1] || observation_oos[cnt] > predict_interval[2]) {
      coverage_binary[cnt] <- 0
    } else {
      coverage_binary[cnt] <- 1
    }
  }
  
  return(mean(coverage_binary))
}


#' Compute empirical coverage rates for individual quantiles
#' 
#' @param observation_oos Out-of-sample observations (vector)
#' @param predict_quantile_mat Matrix of predicted quantiles (rows: observations, cols: quantiles)
#' @param coverage_quantile_vec Vector of quantile levels
#' @return Vector of coverage rates (proportion of observations exceeding each quantile)
compute_all_coverage_rates <- function(observation_oos, predict_quantile_mat, coverage_quantile_vec) {
  coverage_all <- apply(predict_quantile_mat, 2, function(col) {
    mean(observation_oos <= col)
  })
  names(coverage_all) <- coverage_quantile_vec
  return(coverage_all)
}


#' Plot coverage curves against quantile levels
#' 
#' @param coverage_rate Matrix of coverage rates (rows: curves, cols: quantile levels)
#' @param coverage_quantile_vec Vector of quantile levels (x-coordinates)
#' @return None (displays ggplot directly)
plot_coverage_curves <- function(coverage_rate, coverage_quantile_vec, mean_cv_coverage = FALSE, na.rm = FALSE) {
  require(ggplot2)
  require(scales)
  
  # Mean of coverage across folds
  if (mean_cv_coverage == TRUE) {
    coverage_rate = matrix(colMeans(coverage_rate, na.rm = na.rm), nrow=1)
  }
  
  # Ensure inputs are sorted and compatible
  coverage_quantile_vec <- sort(coverage_quantile_vec)
  if (ncol(coverage_rate) != length(coverage_quantile_vec)) {
    stop("Number of columns in coverage_rate must match length of coverage_quantile_vec")
  }
  
  # Adapt data for ggplot
  num_curves <- nrow(coverage_rate)
  num_points <- length(coverage_quantile_vec)
  
  plot_data <- data.frame(
    x = rep(coverage_quantile_vec, num_curves),
    y = as.vector(t(coverage_rate)), # Transpose and flatten the matrix
    group = rep(1:num_curves, each = num_points)
  )
  
  # Create connection data for segments to (0,0) and (1,1)
  connection_data <- data.frame(
    group = factor(1:num_curves),
    x_start = rep(0, num_curves),
    y_start = rep(0, num_curves),
    x_end_start = rep(coverage_quantile_vec[1], num_curves),
    y_end_start = coverage_rate[, 1],
    x_end_end = rep(tail(coverage_quantile_vec, 1), num_curves),
    y_end_end = coverage_rate[, num_points]
  )
  
  # Create the plot
  p <- ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black") + # Diagonal reference line
    geom_point(data = plot_data, aes(x = x, y = y, color = factor(group))) +
    geom_line(data = plot_data, aes(x = x, y = y, group = factor(group), color = factor(group))) + 
    # , stat = "smooth", method = "loess", se = FALSE
    geom_segment(data = connection_data, aes(x = x_start, y = y_start, xend = x_end_start, 
                                             yend = y_end_start, color = group)) +
    geom_segment(data = connection_data, aes(x = 1, y = 1, xend = x_end_end, 
                                             yend = y_end_end, color = group)) +
    xlim(0, 1) +
    ylim(0, 1) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    ggtitle("Coverage Curves vs. Quantile Levels") +
    xlab("Quantile Level") +
    ylab("Empirical Coverage") +
    coord_fixed(ratio = 1) +
    labs(color = "Fold Number")
  
  # Remove legend if there is only one fold
  if (nrow(coverage_rate) == 1) {
    p <- p + theme(legend.position = "none")
  }
  
  print(p)
}


#' Generate 3D plots for quantile surfaces
#' 
#' @param data Training data (matrix)
#' @param observation Training observations (vector)
#' @param predict_func Prediction function
#' @param x_feature Index of feature for x-axis (default: 1)
#' @param y_feature Index of feature for y-axis (default: 2)
#' @param fix_values Vector of values to fix other features (default: means of each feature)
#' @return None (displays plot directly)
generate_3d_plots <- function(data, 
                              observation, 
                              predict_func, 
                              x_feature = 1, 
                              y_feature = 2, 
                              fix_values = NULL,
                              confidence_vec_plot = c(0.05, 0.5, 0.95),
                              xlim_3d = NULL,
                              ylim_3d = NULL,
                              zlim_3d = NULL,
                              alpha = 0.3) {
  require(rgl)
  require(viridis)
  
  # Set default xlim_3d, ylim_3d, and zlim_3d if they are NULL
  if (is.null(xlim_3d)) {
    xlim_3d <- range(data[, x_feature])
  }
  if (is.null(ylim_3d)) {
    ylim_3d <- range(data[, y_feature])
  }
  if (is.null(zlim_3d)) {
    zlim_3d <- range(observation)
  }
  
  if (x_feature < 1 || x_feature > ncol(data) || y_feature < 1 || y_feature > ncol(data)) {
    stop("Invalid feature indices: x_feature and y_feature must be between 1 and ncol(data).")
  }
  if (x_feature == y_feature) {
    stop("x_feature and y_feature must be different.")
  }
  
  if (is.null(fix_values)) {
    fix_values <- apply(data, 2, mean)
  } else if (length(fix_values) != ncol(data)) {
    stop("Length of fix_values must equal ncol(data).")
  }
  
  
  theta_1 <- seq(xlim_3d[1], xlim_3d[2], length.out = 200)
  theta_2 <- seq(ylim_3d[1], ylim_3d[2], length.out = 200)
  z_length <- 200 * 200
  
  z_list <- list()
  for (cnt0 in 1:length(confidence_vec_plot)) {
    z_list[[cnt0]] <- outer(
      theta_1, theta_2,
      FUN = Vectorize(
        function(theta_1, theta_2) {
          input_vec <- fix_values
          input_vec[x_feature] <- theta_1
          input_vec[y_feature] <- theta_2
          predict_func(input_vec, p = confidence_vec_plot[cnt0])
        },
        vectorize.args = c("theta_1", "theta_2")
      )
    )
  }
  
  theta_1_gridsurf <- seq(xlim_3d[1], xlim_3d[2], length.out = 20)
  theta_2_gridsurf <- seq(ylim_3d[1], ylim_3d[2], length.out = 20)
  z_list_gridsurf <- list()
  for (cnt0 in 1:length(confidence_vec_plot)) {
    z_list_gridsurf[[cnt0]] <- outer(
      theta_1_gridsurf, theta_2_gridsurf,
      FUN = Vectorize(
        function(theta_1, theta_2) {
          input_vec <- fix_values
          input_vec[x_feature] <- theta_1
          input_vec[y_feature] <- theta_2
          predict_func(input_vec, p = confidence_vec_plot[cnt0])
        },
        vectorize.args = c("theta_1", "theta_2")
      )
    )
  }
  
  z_combined <- unlist(z_list)
  # jet.colors <- colorRampPalette(c("deepskyblue4", "yellow", "firebrick3"))
  # pal <- jet.colors(100)
  viridis.colors <- colorRampPalette(viridis::viridis(100)) # better visual effect
  pal <- viridis.colors(100) 
  col.ind <- cut(unlist(z_combined), 100)
  
  alpha_3d <- alpha
  for (cnt0 in 1:length(confidence_vec_plot)) {
    if (cnt0 == 1) {
      persp3d(
        theta_1, theta_2, z_list[[cnt0]],
        col = pal[col.ind][(1 + (cnt0 - 1) * z_length):(cnt0 * z_length)],
        xlim = xlim_3d, ylim = ylim_3d, zlim = zlim_3d * 1.3,
        alpha = alpha_3d,
        xlab = colnames(data)[x_feature],
        ylab = colnames(data)[y_feature],
        zlab = "Prediction"
      )
      surface3d(theta_1_gridsurf, theta_2_gridsurf, z_list_gridsurf[[cnt0]], front = 'lines')
    } else {
      persp3d(
        theta_1, theta_2, z_list[[cnt0]],
        col = pal[col.ind][(1 + (cnt0 - 1) * z_length):(cnt0 * z_length)],
        xlim = xlim_3d, ylim = ylim_3d, zlim = zlim_3d * 1.3,
        add = TRUE, alpha = alpha_3d
      )
      surface3d(theta_1_gridsurf, theta_2_gridsurf, z_list_gridsurf[[cnt0]], front = 'lines')
    }
  }
  
  plot3d(
    x = data[, x_feature], y = data[, y_feature], z = observation,
    add = TRUE,
    xlim = xlim_3d, ylim = ylim_3d, zlim = zlim_3d,
    lighting = list(ambient = 0.6, diffuse = 1),
    ltheta = -175, lphi = 70
  )
}

