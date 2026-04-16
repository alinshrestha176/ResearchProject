# ============================================================================
# PRODUCTION SIMULATION - OLS-LASSO FINAL VERSION
# Configuration: 1000 MC iterations × 500 bootstrap samples × 96 cores
# Focus: Comprehensive coverage probability analysis with 5 bootstrap CI methods
# Expected runtime: ~6-7 hours
# MATCHED TO WEE-LASSO FRAMEWORK FOR DIRECT COMPARISON
# ============================================================================

separator <- function(n = 70) paste(rep("=", n), collapse = "")

# ============================================================================
# CONFIGURATION - MATCHED TO WEE-LASSO
# ============================================================================
cat(separator(), "\n")
cat("PRODUCTION MODE: OLS-LASSO FINAL RESULTS\n")
cat("Configuration: 1000-500 with 96 cores\n")
cat(separator(), "\n\n")

N_ITERATIONS <- 1000       # MATCHED to WEE-LASSO
BOOT_ITERATIONS <- 500     # MATCHED to WEE-LASSO
N_LAMBDA <- 15             # MATCHED to WEE-LASSO
C_GRID_SIZE <- 5           # MATCHED to WEE-LASSO
N_CORES <- 96              # MATCHED to WEE-LASSO
MASTER_SEED <- 1762        # MATCHED to WEE-LASSO

cat("Configuration:\n")
cat("  Simulation iterations:", N_ITERATIONS, "\n")
cat("  Bootstrap samples:", BOOT_ITERATIONS, "\n")
cat("  Lambda grid size:", N_LAMBDA, "\n")
cat("  C* grid size:", C_GRID_SIZE, "\n")
cat("  Parallel cores:", N_CORES, "\n")
cat("  Monte Carlo SE: ±", round(sqrt(0.95*0.05/N_ITERATIONS), 4), "\n")
cat("  Expected runtime: ~6-7 hours\n\n")

# ============================================================================
# LOAD LIBRARIES
# ============================================================================
required_packages <- c("MASS", "glmnet", "parallel", "foreach", "doParallel")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "http://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

st <- function(a, b) {
  sign(a) * max(abs(a) - b, 0)
}

scad_threshold <- function(z, lambda, a = 3.7) {
  abs_z <- abs(z)
  sgn <- sign(z)
  
  if (abs_z <= lambda) {
    return(sgn * max(0, abs_z - lambda))
  } else if (abs_z <= a * lambda) {
    return(sgn * ((a - 1) * abs_z - a * lambda) / (a - 2))
  } else {
    return(z)
  }
}

scad_threshold_vec <- Vectorize(scad_threshold, vectorize.args = "z")

# ============================================================================
# BOOTSTRAP CI FUNCTIONS - ALL 5 METHODS (MATCHED TO WEE-LASSO)
# ============================================================================

# 1. Percentile Method
percentile_ci <- function(bootstrap_estimates, point_estimate, alpha = 0.05) {
  ci_lower <- quantile(bootstrap_estimates, alpha/2, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_estimates, 1 - alpha/2, na.rm = TRUE)
  return(list(lower = as.numeric(ci_lower), upper = as.numeric(ci_upper)))
}

# 2. Universal (Reverse Percentile) Method
universal_ci <- function(bootstrap_estimates, point_estimate, alpha = 0.05) {
  deviations <- bootstrap_estimates - point_estimate
  q_lower <- quantile(deviations, alpha/2, na.rm = TRUE)
  q_upper <- quantile(deviations, 1 - alpha/2, na.rm = TRUE)
  
  # Reflect around estimate
  ci_lower <- point_estimate - q_upper
  ci_upper <- point_estimate - q_lower
  
  return(list(lower = as.numeric(ci_lower), upper = as.numeric(ci_upper)))
}

# 3. BCa (Bias-Corrected and Accelerated) Method
bca_ci <- function(bootstrap_estimates, point_estimate, alpha = 0.05) {
  B <- length(bootstrap_estimates)
  
  # Bias-correction factor z0
  prop_less <- sum(bootstrap_estimates < point_estimate, na.rm = TRUE) / B
  prop_less <- max(min(prop_less, 0.9999), 0.0001)
  z0 <- qnorm(prop_less)
  
  # Acceleration factor (simplified - set to 0 for computational efficiency)
  a <- 0
  
  # Adjusted quantile levels
  z_alpha_lower <- qnorm(alpha/2)
  z_alpha_upper <- qnorm(1 - alpha/2)
  
  alpha1 <- pnorm(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
  alpha2 <- pnorm(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))
  
  alpha1 <- max(min(alpha1, 0.9999), 0.0001)
  alpha2 <- max(min(alpha2, 0.9999), 0.0001)
  
  ci_lower <- quantile(bootstrap_estimates, alpha1, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_estimates, alpha2, na.rm = TRUE)
  
  return(list(lower = as.numeric(ci_lower), upper = as.numeric(ci_upper)))
}

# 4. Efron Standard Bias-Corrected Method
efron_standard_ci <- function(bootstrap_estimates, point_estimate, alpha = 0.05) {
  B <- length(bootstrap_estimates)
  
  # Bias-correction factor z0
  prop_less <- sum(bootstrap_estimates < point_estimate, na.rm = TRUE) / B
  prop_less <- max(min(prop_less, 0.9999), 0.0001)
  z0 <- qnorm(prop_less)
  
  # Adjusted quantile levels (no acceleration)
  z_alpha_lower <- qnorm(alpha/2)
  z_alpha_upper <- qnorm(1 - alpha/2)
  
  alpha1 <- pnorm(2*z0 + z_alpha_lower)
  alpha2 <- pnorm(2*z0 + z_alpha_upper)
  
  alpha1 <- max(min(alpha1, 0.9999), 0.0001)
  alpha2 <- max(min(alpha2, 0.9999), 0.0001)
  
  ci_lower <- quantile(bootstrap_estimates, alpha1, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_estimates, alpha2, na.rm = TRUE)
  
  return(list(lower = as.numeric(ci_lower), upper = as.numeric(ci_upper)))
}

# 5. Efron Reflection Method
efron_reflection_ci <- function(bootstrap_estimates, point_estimate, alpha = 0.05) {
  boot_mean <- mean(bootstrap_estimates, na.rm = TRUE)
  q_lower <- quantile(bootstrap_estimates, alpha/2, na.rm = TRUE)
  q_upper <- quantile(bootstrap_estimates, 1 - alpha/2, na.rm = TRUE)
  
  ci_lower <- point_estimate + (boot_mean - q_upper)
  ci_upper <- point_estimate + (boot_mean - q_lower)
  
  return(list(lower = as.numeric(ci_lower), upper = as.numeric(ci_upper)))
}

# ============================================================================
# CORE ESTIMATION FUNCTIONS
# ============================================================================

alpha_estimation <- function(X_list, y_lst, rweights, c_star) {
  X_list <- as.matrix(X_list)
  y_lst <- as.numeric(y_lst)
  rweights <- as.numeric(rweights)
  
  datax <- cbind(1, X_list)
  datax <- as.matrix(datax)
  datay <- matrix(y_lst, ncol = 1)
  
  n <- nrow(datax)
  pdim <- ncol(datax)
  p <- pdim - 1
  
  # Index set determination (first 4 level terms only)
  max_x <- apply(abs(datax[, 2:5]), 2, max)
  threshold <- (n^(-1/2)) * log(n) * max_x
  I_set <- which(threshold < c_star)
  
  I_set_cols <- I_set + 1
  non_I_set_indices <- setdiff(2:5, I_set_cols)
  
  # Weight matrix construction
  wmatdiag <- array(0, c(pdim, pdim, n))
  
  for (t in 1:n) {
    diag(wmatdiag[, , t]) <- 1
    
    if (length(non_I_set_indices) > 0) {
      sum_squares <- sum(datax[t, non_I_set_indices]^2)
      sum_squares <- max(sum_squares, 1e-8)
      
      for (i in non_I_set_indices) {
        wmatdiag[i, i, t] <- 1 / sqrt(sum_squares + 1)
      }
    }
  }
  
  # Weighted estimation
  X_weighted <- matrix(0, nrow = n, ncol = pdim)
  for (i in 1:n) {
    X_weighted[i, ] <- as.numeric(rweights[i]) * (wmatdiag[, , i] %*% datax[i, ])
  }
  
  X_weighted <- matrix(as.numeric(X_weighted), nrow = n, ncol = pdim)
  
  XtX <- t(X_weighted) %*% X_weighted
  ridge_factor <- 1e-8 * mean(diag(XtX))
  diag(XtX) <- diag(XtX) + ridge_factor
  
  Xty <- t(X_weighted) %*% datay
  
  cond_num <- kappa(XtX)
  
  if (cond_num > 1e10) {
    beta <- MASS::ginv(XtX) %*% Xty
  } else {
    beta <- solve(XtX, Xty)
  }
  
  return(list(
    beta = as.numeric(beta),
    I_set = I_set_cols,
    condition_number = cond_num
  ))
}

rho_only_estimation <- function(datalist, rweights) {
  x_1list <- datalist
  n <- length(x_1list)
  x_1t <- x_1list[2:n]
  x_1tless <- x_1list[1:(n-1)]
  
  datax <- matrix(x_1tless, ncol = 1)
  datay <- matrix(x_1t, ncol = 1)
  
  boot_weight <- rweights[2:n]
  wt_matrx <- diag(boot_weight)
  
  XtWX <- t(datax) %*% wt_matrx %*% datax
  XtWy <- t(datax) %*% wt_matrx %*% datay
  
  if(abs(det(XtWX)) < 1e-10) {
    rho_hat <- qr.solve(XtWX, XtWy)
  } else {
    rho_hat <- solve(XtWX, XtWy)
  }
  
  return(as.vector(rho_hat))
}

adaptive_lasso_wls_cd <- function(X, y, lambda, init_beta = NULL, alpha = 1, 
                                  max_iter = 2000, tol = 1e-6, 
                                  include_intercept = TRUE, indices = NULL, 
                                  boot_weight = NULL, mth = 1, ridge_lambda = 0.1) {
  X <- as.matrix(X)
  y <- as.vector(y)
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(boot_weight)) boot_weight <- rep(1, n)
  
  # Center data
  if (include_intercept) {
    y_mean <- sum(boot_weight * y) / sum(boot_weight)
    y <- y - y_mean
    X_means <- colSums(boot_weight * X) / sum(boot_weight)
    X <- scale(X, center = X_means, scale = FALSE)
  } else {
    y_mean <- 0
    X_means <- rep(0, p)
  }
  
  # Initialize beta with ridge regression
  if (is.null(init_beta)) {
    if (all(boot_weight == 1)) {
      XtX <- crossprod(X)
      Xty <- crossprod(X, y)
      ridge_matrix <- XtX + ridge_lambda * diag(p)
      init_beta <- solve(ridge_matrix, Xty)
    } else {
      W <- diag(boot_weight)
      XtWX <- crossprod(X, W %*% X)
      XtWy <- crossprod(X, W %*% y)
      ridge_matrix <- XtWX + ridge_lambda * diag(p)
      init_beta <- solve(ridge_matrix, XtWy)
    }
  } else if (length(init_beta) != p) {
    stop("init_beta must have length equal to number of predictors")
  }
  
  epsilon <- 1e-4
  pen_weights <- 1 / (abs(init_beta) + epsilon)^alpha
  
  beta <- init_beta
  
  # Coordinate descent loop
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    for (j in 1:p) {
      r_j <- y - X %*% beta + X[, j] * beta[j]
      
      # Apply combined weights (mth = 2 for OLS, no adaptive weights)
      if (mth == 1) {
        # WEE method (not used in OLS-LASSO)
        if (j %in% indices) {
          obs_weights <- 1 / (sqrt((X[,j])^2 + 1))
        } else {
          obs_weights <- rep(1, n)
        }
        combined_weights <- boot_weight * obs_weights
      } else if (mth == 2) {
        # OLS method
        combined_weights <- boot_weight
      }
      
      numerator <- sum(combined_weights * X[, j] * r_j)
      denominator <- sum(combined_weights * X[, j]^2)
      z_j <- numerator / denominator
      
      # Adaptive LASSO soft-thresholding
      penalty_scaled <- lambda * pen_weights[j] / denominator
      beta[j] <- sign(z_j) * max(0, abs(z_j) - penalty_scaled)
    }
    
    if (sum(abs(beta - beta_old)) < tol) break
  }
  
  intercept <- y_mean - sum(X_means * beta)
  
  return(list(
    coefficients = if (include_intercept) c(intercept, beta) else beta,
    intercept = intercept,
    beta = beta,
    iterations = iter,
    converged = iter < max_iter,
    lambda = lambda
  ))
}

# ============================================================================
# CROSS-VALIDATION AND LAMBDA GRID
# ============================================================================

ts_cv_5fold <- function(X_list, y_lst, c_values = seq(2, 15, length.out = 10)) {
  
  n <- nrow(X_list)
  fold_size <- floor((n - 500) * 0.2)
  
  results <- matrix(NA, nrow = length(c_values), ncol = 5)
  rownames(results) <- c_values
  
  for (i in seq_along(c_values)) {
    c_current <- c_values[i]
    
    for (fold in 1:5) {
      train_end <- (fold - 1) * fold_size + 500
      test_start <- train_end + 1
      test_end <- min(test_start + fold_size - 1, n)
      
      train_indices <- 1:train_end
      test_indices <- test_start:test_end
      
      X_train <- X_list[train_indices, ]
      y_train <- y_lst[train_indices]
      X_test <- X_list[test_indices, ]
      y_test <- y_lst[test_indices]
      
      rweights <- rep(1, nrow(X_train))
      
      model_fit <- alpha_estimation(X_train, y_train, rweights, c_star = c_current)
      
      X_test1 <- cbind(1, X_test)
      y_pred <- X_test1 %*% model_fit$beta
      
      mse <- mean((y_test - y_pred)^2)
      results[i, fold] <- mse
    }
  }
  
  avg_mse <- rowMeans(results, na.rm = TRUE)
  optimal_c <- c_values[which.min(avg_mse)]
  
  return(list(
    optimal_c = optimal_c,
    all_mse = avg_mse,
    fold_results = results
  ))
}

generate_lambda_grid <- function(X, y, boot_weight = NULL, pen_weights = NULL,
                                 lambda_min_ratio = 0.001, n_lambda = 50) {
  X <- as.matrix(X)
  y <- as.vector(y)
  n <- nrow(X)
  
  if (is.null(boot_weight)) boot_weight <- rep(1, n)
  
  Wy <- boot_weight * y
  
  if (is.null(pen_weights)) {
    lambda_max <- max(abs(t(X) %*% Wy)) / n
  } else {
    lambda_max <- max(abs(t(X) %*% Wy) / pen_weights)
  }
  
  lambda_min <- lambda_min_ratio * lambda_max
  
  grid <- exp(seq(log(lambda_max), log(lambda_min), length.out = n_lambda))
  return(grid)
}

# ============================================================================
# MAIN SIMULATION FUNCTION - OLS-LASSO
# ============================================================================

simulate_one_iteration <- function(g) {
  
  # SET SEED for reproducibility
  set.seed(MASTER_SEED + g)
  
  if (g %% 10 == 0) {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] Iteration", g, "/", N_ITERATIONS, "\n")
  }
  
  library(MASS)
  library(glmnet)
  
  # Data generation
  r <- 25
  n <- 3000
  
  rho <- 0.5
  Sigma <- diag(r)
  
  samples <- mvrnorm(n = n, mu = rep(0, r), Sigma = Sigma)
  df <- data.frame(matrix(nrow = n, ncol = r))
  
  a <- c(0.4, 0.8, 1, 0.99)
  b <- runif((r-4), min = 0, max = 0.7)
  rho <- c(a, b)
  
  for(i in 1:r) {
    df[,i] <- filter(samples[,i], filter = c(rho[i]), method = "recursive", sides = 1, init = 0)
  }
  
  beta <- c(2, 0.5, 1, 1.5, -1)
  error <- rnorm(n)
  
  error1 <- numeric(n)
  error2 <- numeric(n)
  error3 <- numeric(n)
  error4 <- numeric(n)
  
  for (i in 2:n) {
    error1[i] <- df[i, 1] - rho[1] * df[i-1, 1]
    error2[i] <- df[i, 2] - rho[2] * df[i-1, 2]
    error3[i] <- df[i, 3] - rho[3] * df[i-1, 3]
    error4[i] <- df[i, 4] - rho[4] * df[i-1, 4]
  }
  
  y_var <- beta[1] + beta[2] * df[1:(n-1), 1] + beta[3] * df[1:(n-1), 2] +
    beta[4] * df[1:(n-1), 3] + beta[5] * df[1:(n-1), 4] +
    0 * error1[2:n] + 0.5 * error2[2:n] + 1 * error3[2:n] + 0 * error4[2:n] + error[2:n]
  
  X <- cbind(df)
  y_var <- y_var[51:(n-1)]
  X <- X[51:(n-1),]
  
  diff_df <- data.frame(matrix(nrow = n, ncol = r))
  for(i in 1:r) {
    diff_df[[i]] <- c(0, diff(df[,i]))
  }
  
  X <- cbind(X, diff_df[52:n,])
  
  p <- ncol(X)
  new_names <- c(paste0("x", 1:r), paste0("dx", 1:r))
  colnames(X) <- new_names
  dataframe1 <- cbind(y_var, X)
  
  # OLS estimation
  OLS <- lm(y_var~., data = dataframe1)
  OLS_gamma <- OLS$coefficients
  
  x <- model.matrix(y_var~., data = dataframe1)[,-1]
  y <- dataframe1$y_var
  
  X <- cbind(1, x)
  X_with_intercept <- X
  p <- ncol(x)
  n <- nrow(x)
  
  # Cross-validation for c*
  c_grid <- seq(0, 3, length.out = C_GRID_SIZE)
  cv_results <- ts_cv_5fold(X_list = x, y_lst = y, c_values = c_grid)
  optimal_c_star <- cv_results$optimal_c
  
  max_x <- apply(x, 2, function(y) max(abs(y)))
  std <- (n^(-1/2)) * log(n) * max_x
  I_set <- which((n^(-1/2)) * log(n) * max_x < optimal_c_star)
  non_I_set_indices <- setdiff(1:4, I_set)
  
  # Lambda grid for LASSO
  LASSO_grid <- generate_lambda_grid(x, y, n_lambda = N_LAMBDA)
  
  # =========================================================================
  # OLS-LASSO MAIN ESTIMATION
  # =========================================================================
  OLS_LASSO_SIC_values <- numeric(length(LASSO_grid))
  for (j in 1:length(LASSO_grid)) {
    gamma <- adaptive_lasso_wls_cd(x, y, initial_betas = NULL, 
                                   lambda = LASSO_grid[j], 
                                   mth = 2, 
                                   indices = non_I_set_indices)
    LASSO_y_pred <- X %*% gamma$coefficients
    LASSO_RSS <- sum((y - LASSO_y_pred)^2)
    c_n <- max(log(n), log(log(n)) * log(p))
    k <- sum(gamma$beta != 0)
    # FIXED: Use consistent SIC formula (Form 2)
    OLS_LASSO_SIC_values[j] <- n * log(LASSO_RSS/n) + c_n * k
  }
  
  OLS_LASSO_optimal_lambda <- LASSO_grid[which.min(OLS_LASSO_SIC_values)]
  OLS_LASSO_gamma <- adaptive_lasso_wls_cd(x, y, initial_betas = NULL, 
                                           lambda = OLS_LASSO_optimal_lambda, 
                                           mth = 2, 
                                           indices = non_I_set_indices)
  
  # Estimate rho with OLS method
  rweight1 <- rep(1, n)
  OLS_rho <- c()
  for (i in 1:r) {
    rho_OLS.mod <- lm(x[2:n,i]~x[1:(n-1),i]-1)
    rho_OLS_est <- rho_OLS.mod$coef
    OLS_rho <- c(OLS_rho, rho_OLS_est)
  }
  
  OLS_LASSO_alpha <- OLS_LASSO_gamma$coefficients[1:(r+1)]
  OLS_LASSO_eta <- OLS_LASSO_gamma$coefficients[(r+2):(2*r+1)]
  
  OLS_LASSO_rho <- c(0, OLS_rho)
  OLS_LASSO_eta <- c(0, as.array(OLS_LASSO_eta))
  OLS_LASSO_alpha <- as.array(OLS_LASSO_alpha)
  
  OLS_LASSO_beta <- OLS_LASSO_alpha + OLS_LASSO_eta * (OLS_LASSO_rho - 1)
  
  # =========================================================================
  # BOOTSTRAP SECTION
  # =========================================================================
  bootstrap_OLS_LASSO_beta <- matrix(NA, nrow = BOOT_ITERATIONS, ncol = 5)
  bootstrap_OLS_LASSO_eta <- matrix(NA, nrow = BOOT_ITERATIONS, ncol = 4)
  
  for(h in 1:BOOT_ITERATIONS) {
    
    # SET SEED for reproducible bootstrap
    set.seed(MASTER_SEED * 10000 + g * 1000 + h)
    
    boot_weight <- rnorm(n, mean = 1, sd = 1)
    
    # Bootstrap OLS for initial estimates
    W <- diag(boot_weight)
    XtWX <- t(X_with_intercept) %*% W %*% X_with_intercept
    XtWy <- t(X_with_intercept) %*% W %*% y
    
    # FIXED: Stronger ridge regularization for bootstrap stability
    ridge_factor <- 1e-6 * mean(diag(XtWX))
    diag(XtWX) <- diag(XtWX) + ridge_factor
    
    boot_OLS_gamma <- solve(XtWX, XtWy)
    
    # Bootstrap lambda grid
    boot_LASSO_grid <- generate_lambda_grid(x, y, boot_weight = boot_weight, 
                                            n_lambda = N_LAMBDA)
    
    # Bootstrap rho estimation (OLS method)
    boot_OLS_rho <- c()
    for (i in 1:r) {
      rho_OLS_est <- rho_only_estimation(x[,i], boot_weight)
      boot_OLS_rho <- c(boot_OLS_rho, rho_OLS_est)
    }
    boot_OLS_rho <- c(0, boot_OLS_rho)
    
    # OLS-LASSO with bootstrap weights
    OLS_LASSO_SIC_values <- numeric(length(boot_LASSO_grid))
    
    for (j in 1:length(boot_LASSO_grid)) {
      gamma <- adaptive_lasso_wls_cd(x, y, lambda = boot_LASSO_grid[j], 
                                     init_beta = boot_OLS_gamma[-1], 
                                     boot_weight = boot_weight, 
                                     mth = 2, 
                                     indices = non_I_set_indices)
      
      LASSO_y_pred <- X_with_intercept %*% gamma$coefficients
      residuals <- y - LASSO_y_pred
      
      # CRITICAL FIX: Proper weighted RSS normalization
      sum_weights <- sum(boot_weight)
      LASSO_RSS_weighted <- sum(boot_weight * residuals^2) / sum_weights
      
      c_n <- max(log(n), log(log(n)) * log(p))
      k <- sum(gamma$beta != 0)
      # FIXED: Use consistent SIC formula (Form 2) with properly normalized RSS
      OLS_LASSO_SIC_values[j] <- n * log(LASSO_RSS_weighted) + c_n * k
    }
    
    OLS_LASSO_optimal_lambda <- boot_LASSO_grid[which.min(OLS_LASSO_SIC_values)]
    
    boot_OLS_LASSO_gamma <- adaptive_lasso_wls_cd(x, y, 
                                                  lambda = OLS_LASSO_optimal_lambda, 
                                                  init_beta = boot_OLS_gamma[-1], 
                                                  boot_weight = boot_weight, 
                                                  mth = 2, 
                                                  indices = non_I_set_indices)
    
    boot_OLS_LASSO_alpha <- boot_OLS_LASSO_gamma$coefficients[1:(r+1)]
    boot_OLS_LASSO_eta <- boot_OLS_LASSO_gamma$coefficients[(r+2):(2*r+1)]
    boot_OLS_LASSO_eta <- c(0, as.array(boot_OLS_LASSO_eta))
    boot_OLS_LASSO_alpha <- as.array(boot_OLS_LASSO_alpha)
    
    boot_OLS_LASSO_beta <- boot_OLS_LASSO_alpha + boot_OLS_LASSO_eta * (boot_OLS_rho - 1)
    
    bootstrap_OLS_LASSO_beta[h,] <- boot_OLS_LASSO_beta[1:5]
    bootstrap_OLS_LASSO_eta[h,] <- boot_OLS_LASSO_eta[2:5]
  }
  
  # =========================================================================
  # COMPUTE ALL 5 CI METHODS FOR BETA
  # =========================================================================
  
  # 1. Percentile Method
  OLS_LASSO_beta_CI_percentile <- matrix(NA, nrow = 2, ncol = 5)
  for (idx in 1:5) {
    ci <- percentile_ci(bootstrap_OLS_LASSO_beta[, idx], OLS_LASSO_beta[idx])
    OLS_LASSO_beta_CI_percentile[1, idx] <- ci$lower
    OLS_LASSO_beta_CI_percentile[2, idx] <- ci$upper
  }
  
  # 2. Universal Method
  OLS_LASSO_beta_CI_universal <- matrix(NA, nrow = 2, ncol = 5)
  for (idx in 1:5) {
    ci <- universal_ci(bootstrap_OLS_LASSO_beta[, idx], OLS_LASSO_beta[idx])
    OLS_LASSO_beta_CI_universal[1, idx] <- ci$lower
    OLS_LASSO_beta_CI_universal[2, idx] <- ci$upper
  }
  
  # 3. BCa Method
  OLS_LASSO_beta_CI_bca <- matrix(NA, nrow = 2, ncol = 5)
  for (idx in 1:5) {
    ci <- bca_ci(bootstrap_OLS_LASSO_beta[, idx], OLS_LASSO_beta[idx])
    OLS_LASSO_beta_CI_bca[1, idx] <- ci$lower
    OLS_LASSO_beta_CI_bca[2, idx] <- ci$upper
  }
  
  # 4. Efron Standard Method
  OLS_LASSO_beta_CI_efron_standard <- matrix(NA, nrow = 2, ncol = 5)
  for (idx in 1:5) {
    ci <- efron_standard_ci(bootstrap_OLS_LASSO_beta[, idx], OLS_LASSO_beta[idx])
    OLS_LASSO_beta_CI_efron_standard[1, idx] <- ci$lower
    OLS_LASSO_beta_CI_efron_standard[2, idx] <- ci$upper
  }
  
  # 5. Efron Reflection Method
  OLS_LASSO_beta_CI_efron_reflection <- matrix(NA, nrow = 2, ncol = 5)
  for (idx in 1:5) {
    ci <- efron_reflection_ci(bootstrap_OLS_LASSO_beta[, idx], OLS_LASSO_beta[idx])
    OLS_LASSO_beta_CI_efron_reflection[1, idx] <- ci$lower
    OLS_LASSO_beta_CI_efron_reflection[2, idx] <- ci$upper
  }
  
  # For ETA: use percentile method
  OLS_LASSO_eta_CI_percentile <- matrix(NA, nrow = 2, ncol = 4)
  for (idx in 1:4) {
    ci <- percentile_ci(bootstrap_OLS_LASSO_eta[, idx], OLS_LASSO_eta[idx+1])
    OLS_LASSO_eta_CI_percentile[1, idx] <- ci$lower
    OLS_LASSO_eta_CI_percentile[2, idx] <- ci$upper
  }
  
  # =========================================================================
  # COMPUTE COVERAGE FOR ALL METHODS
  # =========================================================================
  
  true_beta_vec <- c(2, 0.5, 1, 1.5, -1)
  true_eta_vec <- c(0, 0.5, 1, 0)
  
  OLS_LASSO_beta_coverage_percentile <- (OLS_LASSO_beta_CI_percentile[1,] <= true_beta_vec) & 
    (OLS_LASSO_beta_CI_percentile[2,] >= true_beta_vec)
  
  OLS_LASSO_beta_coverage_universal <- (OLS_LASSO_beta_CI_universal[1,] <= true_beta_vec) & 
    (OLS_LASSO_beta_CI_universal[2,] >= true_beta_vec)
  
  OLS_LASSO_beta_coverage_bca <- (OLS_LASSO_beta_CI_bca[1,] <= true_beta_vec) & 
    (OLS_LASSO_beta_CI_bca[2,] >= true_beta_vec)
  
  OLS_LASSO_beta_coverage_efron_standard <- (OLS_LASSO_beta_CI_efron_standard[1,] <= true_beta_vec) & 
    (OLS_LASSO_beta_CI_efron_standard[2,] >= true_beta_vec)
  
  OLS_LASSO_beta_coverage_efron_reflection <- (OLS_LASSO_beta_CI_efron_reflection[1,] <= true_beta_vec) & 
    (OLS_LASSO_beta_CI_efron_reflection[2,] >= true_beta_vec)
  
  OLS_LASSO_eta_coverage <- (OLS_LASSO_eta_CI_percentile[1,] <= true_eta_vec) & 
    (OLS_LASSO_eta_CI_percentile[2,] >= true_eta_vec)
  
  # =========================================================================
  # VARIABLE SELECTION METRICS
  # =========================================================================
  
  true_beta_nonzero <- true_beta_vec != 0
  true_eta_nonzero <- true_eta_vec != 0
  
  OLS_LASSO_beta_fn <- sum(true_beta_nonzero & (abs(OLS_LASSO_beta[1:5]) < 1e-6))
  OLS_LASSO_beta_fp <- sum(!true_beta_nonzero & (abs(OLS_LASSO_beta[1:5]) >= 1e-6))
  OLS_LASSO_eta_fn <- sum(true_eta_nonzero & (abs(OLS_LASSO_eta[2:5]) < 1e-6))
  OLS_LASSO_eta_fp <- sum(!true_eta_nonzero & (abs(OLS_LASSO_eta[2:5]) >= 1e-6))
  
  # =========================================================================
  # BOOTSTRAP STATISTICS
  # =========================================================================
  
  OLS_LASSO_beta_boot_mean <- colMeans(bootstrap_OLS_LASSO_beta, na.rm = TRUE)
  OLS_LASSO_beta_boot_sd <- apply(bootstrap_OLS_LASSO_beta, 2, sd, na.rm = TRUE)
  OLS_LASSO_eta_boot_mean <- colMeans(bootstrap_OLS_LASSO_eta, na.rm = TRUE)
  OLS_LASSO_eta_boot_sd <- apply(bootstrap_OLS_LASSO_eta, 2, sd, na.rm = TRUE)
  
  return(list(
    OLS_LASSO_beta = OLS_LASSO_beta[1:5],
    OLS_LASSO_eta = OLS_LASSO_eta[2:5],
    
    # Coverage for all five methods
    OLS_LASSO_beta_coverage_percentile = OLS_LASSO_beta_coverage_percentile,
    OLS_LASSO_beta_coverage_universal = OLS_LASSO_beta_coverage_universal,
    OLS_LASSO_beta_coverage_bca = OLS_LASSO_beta_coverage_bca,
    OLS_LASSO_beta_coverage_efron_standard = OLS_LASSO_beta_coverage_efron_standard,
    OLS_LASSO_beta_coverage_efron_reflection = OLS_LASSO_beta_coverage_efron_reflection,
    OLS_LASSO_eta_coverage = OLS_LASSO_eta_coverage,
    
    # CI bounds for all methods
    OLS_LASSO_beta_CI_percentile_lower = OLS_LASSO_beta_CI_percentile[1, ],
    OLS_LASSO_beta_CI_percentile_upper = OLS_LASSO_beta_CI_percentile[2, ],
    OLS_LASSO_beta_CI_universal_lower = OLS_LASSO_beta_CI_universal[1, ],
    OLS_LASSO_beta_CI_universal_upper = OLS_LASSO_beta_CI_universal[2, ],
    OLS_LASSO_beta_CI_bca_lower = OLS_LASSO_beta_CI_bca[1, ],
    OLS_LASSO_beta_CI_bca_upper = OLS_LASSO_beta_CI_bca[2, ],
    OLS_LASSO_beta_CI_efron_standard_lower = OLS_LASSO_beta_CI_efron_standard[1, ],
    OLS_LASSO_beta_CI_efron_standard_upper = OLS_LASSO_beta_CI_efron_standard[2, ],
    OLS_LASSO_beta_CI_efron_reflection_lower = OLS_LASSO_beta_CI_efron_reflection[1, ],
    OLS_LASSO_beta_CI_efron_reflection_upper = OLS_LASSO_beta_CI_efron_reflection[2, ],
    
    # Variable selection metrics
    OLS_LASSO_beta_fn = OLS_LASSO_beta_fn,
    OLS_LASSO_beta_fp = OLS_LASSO_beta_fp,
    OLS_LASSO_eta_fn = OLS_LASSO_eta_fn,
    OLS_LASSO_eta_fp = OLS_LASSO_eta_fp,
    
    # Bootstrap statistics
    OLS_LASSO_beta_boot_mean = OLS_LASSO_beta_boot_mean,
    OLS_LASSO_beta_boot_sd = OLS_LASSO_beta_boot_sd,
    OLS_LASSO_eta_boot_mean = OLS_LASSO_eta_boot_mean,
    OLS_LASSO_eta_boot_sd = OLS_LASSO_eta_boot_sd
  ))
  
}

# ============================================================================
# PARALLEL EXECUTION
# ============================================================================

cat("\n", separator(), "\n")
cat("SETTING UP PARALLEL PROCESSING\n")
cat(separator(), "\n\n")

cl <- makeCluster(N_CORES)
registerDoParallel(cl)

clusterExport(cl, c(
  "simulate_one_iteration",
  "alpha_estimation",
  "rho_only_estimation", 
  "adaptive_lasso_wls_cd",
  "ts_cv_5fold",
  "st",
  "scad_threshold",
  "scad_threshold_vec",
  "generate_lambda_grid",
  "percentile_ci",
  "universal_ci",
  "bca_ci",
  "efron_standard_ci",
  "efron_reflection_ci",
  "N_ITERATIONS",
  "BOOT_ITERATIONS",
  "N_LAMBDA",
  "C_GRID_SIZE",
  "MASTER_SEED"
))

clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
})

cat("\n", separator(), "\n")
cat("STARTING SIMULATION\n")
cat(separator(), "\n\n")

start_time <- Sys.time()

results <- foreach(
  i = 1:N_ITERATIONS,
  .packages = c("MASS", "glmnet"),
  .errorhandling = "pass"
) %dopar% {
  simulate_one_iteration(i)
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "hours")

stopCluster(cl)

cat("\n", separator(), "\n")
cat("SIMULATION COMPLETE\n")
cat("Total time:", round(total_time, 2), "hours\n")
cat(separator(), "\n\n")

# ============================================================================
# RESULTS ANALYSIS
# ============================================================================

cat("Aggregating results...\n\n")

# Check for errors
error_indices <- which(sapply(results, function(x) inherits(x, "error")))
if(length(error_indices) > 0) {
  cat("WARNING:", length(error_indices), "iterations failed\n\n")
  results <- results[-error_indices]
}

true_beta <- c(2, 0.5, 1, 1.5, -1)
true_eta <- c(0, 0.5, 1, 0)

# For BETA values
OLS_LASSO_beta_df <- do.call(rbind, lapply(results, function(x) x$OLS_LASSO_beta))
OLS_LASSO_beta_mean <- round(colMeans(OLS_LASSO_beta_df), 3)
OLS_LASSO_beta_sd <- round(apply(OLS_LASSO_beta_df, 2, sd), 3)

# For ETA values
OLS_LASSO_eta_df <- do.call(rbind, lapply(results, function(x) x$OLS_LASSO_eta))
OLS_LASSO_eta_mean <- round(colMeans(OLS_LASSO_eta_df), 3)
OLS_LASSO_eta_sd <- round(apply(OLS_LASSO_eta_df, 2, sd), 3)

cat("OLS_LASSO_beta_mean:", paste(OLS_LASSO_beta_mean, collapse = ", "), "\n")
cat("OLS_LASSO_beta_sd:", paste(OLS_LASSO_beta_sd, collapse = ", "), "\n\n")

cat("OLS_LASSO_eta_mean:", paste(OLS_LASSO_eta_mean, collapse = ", "), "\n")
cat("OLS_LASSO_eta_sd:", paste(OLS_LASSO_eta_sd, collapse = ", "), "\n\n")

# ============================================================================
# COVERAGE PROBABILITIES - ALL 5 METHODS
# ============================================================================

cat("\n", separator(110), "\n")
cat("OLS-LASSO COVERAGE PROBABILITY ANALYSIS (ALL 5 METHODS)\n")
cat(separator(110), "\n\n")

OLS_LASSO_beta_coverage_percentile_prob <- colMeans(do.call(rbind, lapply(results, function(x) x$OLS_LASSO_beta_coverage_percentile)))
OLS_LASSO_beta_coverage_universal_prob <- colMeans(do.call(rbind, lapply(results, function(x) x$OLS_LASSO_beta_coverage_universal)))
OLS_LASSO_beta_coverage_bca_prob <- colMeans(do.call(rbind, lapply(results, function(x) x$OLS_LASSO_beta_coverage_bca)))
OLS_LASSO_beta_coverage_efron_standard_prob <- colMeans(do.call(rbind, lapply(results, function(x) x$OLS_LASSO_beta_coverage_efron_standard)))
OLS_LASSO_beta_coverage_efron_reflection_prob <- colMeans(do.call(rbind, lapply(results, function(x) x$OLS_LASSO_beta_coverage_efron_reflection)))
OLS_LASSO_eta_coverage_prob <- colMeans(do.call(rbind, lapply(results, function(x) x$OLS_LASSO_eta_coverage)))

cat("=== BETA COVERAGE PROBABILITIES (Target: 0.95) ===\n\n")
cat(sprintf("%-10s | %-11s | %-11s | %-11s | %-13s | %-13s\n", 
            "Parameter", "Percentile", "Universal", "BCa", "Efron-Std", "Efron-Refl"))
cat(sprintf("%-10s | %-11s | %-11s | %-11s | %-13s | %-13s\n", 
            "", "(no corr.)", "(reflect)", "(adaptive)", "(quantile)", "(mean shift)"))
cat(paste(rep("-", 110), collapse=""), "\n")

for(i in 1:5) {
  cat(sprintf("Beta[%d]    | %11.3f | %11.3f | %11.3f | %13.3f | %13.3f", 
              i, 
              OLS_LASSO_beta_coverage_percentile_prob[i],
              OLS_LASSO_beta_coverage_universal_prob[i],
              OLS_LASSO_beta_coverage_bca_prob[i],
              OLS_LASSO_beta_coverage_efron_standard_prob[i],
              OLS_LASSO_beta_coverage_efron_reflection_prob[i]))
  
  if (i == 2) {
    cat("  <- FOCUS")
  }
  cat("\n")
}
cat("\n")

cat("=== ETA COVERAGE PROBABILITIES (Percentile Method, Target: 0.95) ===\n\n")
for(i in 1:4) {
  cat("Eta[", i, "]: ", sprintf("%.3f", OLS_LASSO_eta_coverage_prob[i]), "\n", sep="")
}
cat("\n")

# ============================================================================
# METHOD COMPARISON
# ============================================================================

cat(separator(110), "\n")
cat("METHOD COMPARISON SUMMARY\n")
cat(separator(110), "\n\n")

# Count parameters with adequate coverage (>= 0.935)
percentile_pass <- sum(OLS_LASSO_beta_coverage_percentile_prob >= 0.935)
universal_pass <- sum(OLS_LASSO_beta_coverage_universal_prob >= 0.935)
bca_pass <- sum(OLS_LASSO_beta_coverage_bca_prob >= 0.935)
efron_std_pass <- sum(OLS_LASSO_beta_coverage_efron_standard_prob >= 0.935)
efron_refl_pass <- sum(OLS_LASSO_beta_coverage_efron_reflection_prob >= 0.935)

cat("Parameters with adequate coverage (>=0.935):\n")
cat("  Percentile method:        ", percentile_pass, "/ 5\n")
cat("  Universal method:         ", universal_pass, "/ 5\n")
cat("  BCa method:               ", bca_pass, "/ 5\n")
cat("  Efron Standard method:    ", efron_std_pass, "/ 5\n")
cat("  Efron Reflection method:  ", efron_refl_pass, "/ 5\n\n")

# Average absolute deviation from 0.95
percentile_dev <- mean(abs(OLS_LASSO_beta_coverage_percentile_prob - 0.95))
universal_dev <- mean(abs(OLS_LASSO_beta_coverage_universal_prob - 0.95))
bca_dev <- mean(abs(OLS_LASSO_beta_coverage_bca_prob - 0.95))
efron_std_dev <- mean(abs(OLS_LASSO_beta_coverage_efron_standard_prob - 0.95))
efron_refl_dev <- mean(abs(OLS_LASSO_beta_coverage_efron_reflection_prob - 0.95))

cat("Average absolute deviation from 0.95:\n")
cat("  Percentile method:        ", sprintf("%.4f", percentile_dev), "\n")
cat("  Universal method:         ", sprintf("%.4f", universal_dev), "\n")
cat("  BCa method:               ", sprintf("%.4f", bca_dev), "\n")
cat("  Efron Standard method:    ", sprintf("%.4f", efron_std_dev), "\n")
cat("  Efron Reflection method:  ", sprintf("%.4f", efron_refl_dev), "\n\n")

# Best method
all_devs <- c(Percentile = percentile_dev, 
              Universal = universal_dev,
              BCa = bca_dev, 
              `Efron-Std` = efron_std_dev, 
              `Efron-Refl` = efron_refl_dev)
best_method <- names(which.min(all_devs))

cat("Best performing method overall: ", best_method, "\n\n")

# ============================================================================
# BOOTSTRAP DIAGNOSTICS
# ============================================================================

cat(separator(110), "\n")
cat("BOOTSTRAP DIAGNOSTICS\n")
cat(separator(110), "\n\n")

# Average bootstrap mean and SD across simulations
OLS_LASSO_beta_boot_mean_avg <- rowMeans(do.call(cbind, lapply(results, function(x) x$OLS_LASSO_beta_boot_mean)))
OLS_LASSO_beta_boot_sd_avg <- rowMeans(do.call(cbind, lapply(results, function(x) x$OLS_LASSO_beta_boot_sd)))

# Calculate bias
beta_bias <- OLS_LASSO_beta_boot_mean_avg - true_beta

cat("=== BIAS ANALYSIS FOR BETA PARAMETERS ===\n\n")
cat(sprintf("%-10s | %-10s | %-10s | %-10s | %-10s\n", "Parameter", "True", "Boot Mean", "Bias", "Bias/SD"))
cat(paste(rep("-", 60), collapse=""), "\n")
for(i in 1:5) {
  bias_ratio <- beta_bias[i] / OLS_LASSO_beta_boot_sd_avg[i]
  cat(sprintf("Beta[%d]    | %10.3f | %10.3f | %10.4f | %10.4f\n", 
              i, true_beta[i], OLS_LASSO_beta_boot_mean_avg[i], 
              beta_bias[i], bias_ratio))
}
cat("\n")

# Report parameters with poor coverage
cat("=== PARAMETERS WITH POOR COVERAGE (<0.935) IN ANY METHOD ===\n\n")

for(i in 1:5) {
  coverages <- c(
    OLS_LASSO_beta_coverage_percentile_prob[i],
    OLS_LASSO_beta_coverage_universal_prob[i],
    OLS_LASSO_beta_coverage_bca_prob[i],
    OLS_LASSO_beta_coverage_efron_standard_prob[i],
    OLS_LASSO_beta_coverage_efron_reflection_prob[i]
  )
  
  if(any(coverages < 0.935)) {
    cat("   Beta[", i, "]:\n", sep="")
    cat("     Percentile:        ", sprintf("%.3f", coverages[1]), 
        if(coverages[1] < 0.935) " FAIL" else " PASS", "\n")
    cat("     Universal:         ", sprintf("%.3f", coverages[2]), 
        if(coverages[2] < 0.935) " FAIL" else " PASS", "\n")
    cat("     BCa:               ", sprintf("%.3f", coverages[3]), 
        if(coverages[3] < 0.935) " FAIL" else " PASS", "\n")
    cat("     Efron Standard:    ", sprintf("%.3f", coverages[4]), 
        if(coverages[4] < 0.935) " FAIL" else " PASS", "\n")
    cat("     Efron Reflection:  ", sprintf("%.3f", coverages[5]), 
        if(coverages[5] < 0.935) " FAIL" else " PASS", "\n")
    cat("     True value:        ", sprintf("%.3f", true_beta[i]), "\n")
    cat("     Boot mean:         ", sprintf("%.3f", OLS_LASSO_beta_boot_mean_avg[i]), "\n")
    cat("     Boot SD:           ", sprintf("%.3f", OLS_LASSO_beta_boot_sd_avg[i]), "\n")
    cat("     Bias:              ", sprintf("%.4f", beta_bias[i]), "\n")
    cat("     Bias/SD:           ", sprintf("%.4f", beta_bias[i]/OLS_LASSO_beta_boot_sd_avg[i]), "\n\n")
  }
}

# ETA parameters with poor coverage
poor_eta_indices <- which(OLS_LASSO_eta_coverage_prob < 0.935)
if(length(poor_eta_indices) > 0) {
  cat("ETA parameters with poor coverage:\n\n")
  OLS_LASSO_eta_boot_mean_avg <- rowMeans(do.call(cbind, lapply(results, function(x) x$OLS_LASSO_eta_boot_mean)))
  OLS_LASSO_eta_boot_sd_avg <- rowMeans(do.call(cbind, lapply(results, function(x) x$OLS_LASSO_eta_boot_sd)))
  
  for(idx in poor_eta_indices) {
    cat("   Eta[", idx, "]:\n", sep="")
    cat("     Coverage:    ", sprintf("%.3f", OLS_LASSO_eta_coverage_prob[idx]), "\n")
    cat("     True value:  ", sprintf("%.3f", true_eta[idx]), "\n")
    cat("     Boot mean:   ", sprintf("%.3f", OLS_LASSO_eta_boot_mean_avg[idx]), "\n")
    cat("     Boot SD:     ", sprintf("%.3f", OLS_LASSO_eta_boot_sd_avg[idx]), "\n\n")
  }
} else {
  cat("✓ All ETA parameters have adequate coverage (≥0.935)\n\n")
}

# ============================================================================
# ADDITIONAL STATISTICS
# ============================================================================

cat(separator(110), "\n")
cat("ADDITIONAL STATISTICS\n")
cat(separator(110), "\n\n")

# MSE
mse <- function(x, y_true) mean((x - y_true)^2)
OLS_LASSO_beta_mse <- round(mapply(mse, as.data.frame(OLS_LASSO_beta_df), true_beta), 4)

cat("BETA MSE:\n")
for(i in 1:5) {
  cat("Beta[", i, "]: ", OLS_LASSO_beta_mse[i], "\n", sep="")
}
cat("\n")

# Variable selection
OLS_LASSO_beta_FN_total <- sum(sapply(results, function(x) x$OLS_LASSO_beta_fn))
OLS_LASSO_beta_FP_total <- sum(sapply(results, function(x) x$OLS_LASSO_beta_fp))
OLS_LASSO_eta_FN_total <- sum(sapply(results, function(x) x$OLS_LASSO_eta_fn))
OLS_LASSO_eta_FP_total <- sum(sapply(results, function(x) x$OLS_LASSO_eta_fp))

cat("VARIABLE SELECTION:\n")
cat("BETA - False Negatives:", OLS_LASSO_beta_FN_total, ", False Positives:", OLS_LASSO_beta_FP_total, "\n")
cat("ETA  - False Negatives:", OLS_LASSO_eta_FN_total, ", False Positives:", OLS_LASSO_eta_FP_total, "\n\n")

cat(separator(), "\n")
cat("ANALYSIS COMPLETE\n")
cat(separator(), "\n")
