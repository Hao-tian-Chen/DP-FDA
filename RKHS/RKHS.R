############################## RKHS Gaussian Kernel Case #################################

library("doParallel")
library("foreach")

######## Part 1. Data Preparation and Helpful Functions

#### Define the RBF (Gaussian) Kernel Function (Matrix) ####
rbf_kernel <- function(x, y, sigma) {
  # Computes the pairwise RBF kernel matrix between vectors x and y
  dist_sq <- outer(x, y, "-")^2
  exp(-dist_sq / (2 * sigma^2))
}

#### Define the B' choices ####
get_Bprime <- function(n, B, choice) {
  switch(
    choice,
    B = B,
    sqrtlog = sqrt(log(n)),
    log = log(n),
    sqrtn = sqrt(n),
    n_2_3 = n^(2/3),
    stop("Unknown Bprime choice.")
  )
}

#### Latent RKHS coefficient generation ####
generate_rkhs_latent <- function(
    n, M, B, sigma_kernel, mu0 = NULL
) {
  
  c_grid <- seq(0, 1, length.out = M)
  
  G <- rbf_kernel(c_grid, c_grid, sigma_kernel)
  G <- (G + t(G)) / 2
  
  # Small numerical ridge for stable inverse square root
  G_reg <- G + diag(1e-8, M)
  
  eig <- eigen(G_reg, symmetric = TRUE)
  evals <- pmax(eig$values, 1e-12)
  
  G_inv_sqrt <- eig$vectors %*%
    diag(1 / sqrt(evals), nrow = M) %*%
    t(eig$vectors)
  
  # Latent Gaussian mean; distinct from privacy parameter mu
  if (is.null(mu0)) {
    mu0 <- seq(-M / 2, M / 2, length.out = M)
  }
  
  if (length(mu0) != M) {
    stop("mu0 must have length M.")
  }
  
  Z <- matrix(rnorm(n * M), nrow = n, ncol = M)
  Z <- sweep(Z, 2, mu0, "+")
  
  Z_norm <- sqrt(rowSums(Z^2))
  
  # Guard against a numerically zero row norm
  if (any(Z_norm <= .Machine$double.eps)) {
    stop("A latent Gaussian vector has numerically zero norm.")
  }
  
  U <- Z / Z_norm
  
  # Uniform volume distribution inside the M-dimensional ball
  r <- B * runif(n)^(1 / M)
  V <- U * r
  
  # RKHS coefficients: alpha_i' G alpha_i = ||V_i||^2
  alpha_mat <- V %*% t(G_inv_sqrt)
  
  list(
    c_grid = c_grid,
    G = G,
    alpha_mat = alpha_mat,
    alpha_bar_true = colMeans(alpha_mat)
  )
}

#### Evaluate fixed latent RKHS functions on a time grid ####
evaluate_rkhs_data <- function(
    alpha_mat, T_pts, sigma_kernel
) {
  
  t_grid <- seq(0, 1, length.out = T_pts)
  
  M <- ncol(alpha_mat)
  c_grid <- seq(0, 1, length.out = M)
  
  K_TC <- rbf_kernel(
    t_grid, c_grid, sigma_kernel
  )
  
  f_mat <- alpha_mat %*% t(K_TC)
  
  list(
    t_grid = t_grid,
    c_grid = c_grid,
    K_TC = K_TC,
    f_mat = f_mat,
    f_bar = colMeans(f_mat)
  )
}

#### Add observation noise; latent functions stay fixed ####
generate_noisy_dataset <- function(f_mat, sigma_0) {
  
  n <- nrow(f_mat)
  T_pts <- ncol(f_mat)
  
  noise_mat <- matrix(rnorm(n * T_pts, mean = 0, sd = sigma_0),
                      nrow = n, ncol = T_pts)

  list(
    x_mat = f_mat + noise_mat
  )
}

#### Fixed kernel design and eigendecomposition ####
prepare_rkhs_design <- function(
    T_pts, sigma_kernel
) {
  
  t_grid <- seq(0, 1, length.out = T_pts)
  
  K_T <- rbf_kernel(
    t_grid, t_grid, sigma_kernel
  )
  
  K_T <- (K_T + t(K_T)) / 2
  
  eig <- eigen(
    K_T,
    symmetric = TRUE
  )
  
  evals <- pmax(eig$values, 0)
  evecs <- eig$vectors
  
  # For GP DP noise, the RBF matrix should be positive
  # semidefinite; clip tiny negative numerical eigenvalues.
  K_psd <- evecs %*%
    diag(evals, nrow = T_pts) %*%
    t(evecs)
  
  # Cholesky may fail for numerically semidefinite matrices.
  # Use the eigenbasis for robust Gaussian sampling instead.
  
  list(
    t_grid = t_grid,
    K_T = K_T,
    K_psd = K_psd,
    evals = evals,
    evecs = evecs,
    T_pts = T_pts,
    sigma_kernel = sigma_kernel
  )
}


######## Part 2. Mechanism Design 

######## Mechanism A1: FTP with radial clamping ########

run_FTP_regime_clamp <- function(
    x_mat, design, B_prime, n, mu, lambda_min = 1e-8
) {
  
  K_T <- design$K_T
  T_pts <- design$T_pts
  
  # Individual regularized kernel fits
  beta_mat <- t(
    solve(
      K_T + diag(lambda_min, T_pts),
      t(x_mat)
    )
  )
  
  # RKHS norm: sqrt(beta_i' K_T beta_i)
  norms <- sqrt(
    pmax(
      rowSums((beta_mat %*% K_T) * beta_mat),
      0
    )
  )
  
  scale_factors <- rep(1, n)
  idx <- norms > B_prime
  
  scale_factors[idx] <-
    B_prime / norms[idx]
  
  beta_clamped <- beta_mat * scale_factors
  
  beta_mean <- colMeans(beta_clamped)
  
  f_mean_eval <- as.vector(
    K_T %*% beta_mean
  )
  
  # Fresh GP DP noise
  sigma_DP <- 2 * B_prime / (n * mu)
  
  Z <- as.vector(
    design$evecs %*%
      (
        sqrt(design$evals) *
          rnorm(T_pts)
      )
  ) * sigma_DP
  
  list(
    f_noisy_T = f_mean_eval + Z,
    num_clamped = sum(idx),
    prop_clamped = mean(idx),
    norms = norms,
    sigma_DP = sigma_DP
  )
}


######## Mechanism A3: FTP with exact constrained optimization ########

run_FTP_regime_constrained <- function(
    x_mat, design, B_prime, n, mu, lambda_min = 1e-8
) {
  
  K_T <- design$K_T
  evals <- design$evals
  evecs <- design$evecs
  
  T_pts <- design$T_pts
  n_obs <- nrow(x_mat)
  
  # Project observations onto the fixed eigenbasis
  C_mat <- t(evecs) %*% t(x_mat)
  
  beta_mat <- matrix(
    0, nrow = T_pts, ncol = n_obs
  )
  
  unconstrained_norms <- numeric(n_obs)
  lambda_stars <- numeric(n_obs)
  was_constrained <- logical(n_obs)
  
  for (i in seq_len(n_obs)) {
    
    c_vec <- C_mat[, i]
    
    calc_norm_sq <- function(lam) {
      sum(
        evals * c_vec^2 /
          (evals + lam)^2
      )
    }
    
    norm_sq <- calc_norm_sq(lambda_min)
    
    unconstrained_norms[i] <- sqrt(
      max(norm_sq, 0)
    )
    
    if (norm_sq <= B_prime^2) {
      
      lambda_star <- lambda_min
      
    } else {
      
      norm_func <- function(lam) {
        calc_norm_sq(lam) - B_prime^2
      }
      
      lam_upper <- max(1, lambda_min * 2)
      
      while (norm_func(lam_upper) > 0) {
        lam_upper <- lam_upper * 10
      }
      
      lambda_star <- uniroot(
        norm_func,
        lower = lambda_min,
        upper = lam_upper
      )$root
      
      was_constrained[i] <- TRUE
    }
    
    lambda_stars[i] <- lambda_star
    
    beta_mat[, i] <- evecs %*%
      (
        c_vec / (evals + lambda_star)
      )
  }
  
  beta_mean <- rowMeans(beta_mat)
  
  f_mean_eval <- as.vector(
    K_T %*% beta_mean
  )
  
  sigma_DP <- 2 * B_prime / (n * mu)
  
  Z <- sigma_DP * as.vector(
    design$evecs %*%
      (
        sqrt(design$evals) *
          rnorm(T_pts)
      )
  )
  
  list(
    f_noisy_T = f_mean_eval + Z,
    num_constrained = sum(was_constrained),
    prop_constrained = mean(was_constrained),
    unconstrained_norms = unconstrained_norms,
    lambda_stars = lambda_stars,
    sigma_DP = sigma_DP
  )
}


#### PTF preprocessing and DP randomization ####

prepare_PTF_data <- function(
    x_mat, B_prime, n, mu
) {
  
  x_mat_clipped <- pmin(
    pmax(x_mat, -B_prime),
    B_prime
  )
  
  x_bar <- colMeans(x_mat_clipped)
  
  T_pts <- ncol(x_mat)
  
  sigma_j <- 2 * B_prime * sqrt(T_pts) /
    (n * mu)
  
  x_tilde <- x_bar + rnorm(
    T_pts, mean = 0, sd = sigma_j
  )
  
  list(
    x_tilde = x_tilde,
    num_clamped = sum(abs(x_mat) > B_prime),
    prop_clamped = mean(abs(x_mat) > B_prime),
    sigma_j = sigma_j
  )
}

######## Mechanism A2: PTF with constrained KRR ########
fit_PTF_constrained <- function(
    x_tilde, design, B_prime,
    lambda_min = 1e-8
) {
  
  evals <- design$evals
  evecs <- design$evecs
  K_T <- design$K_T
  T_pts <- design$T_pts
  
  u_vec <- as.vector(
    t(evecs) %*% x_tilde
  )
  
  calc_norm_sq <- function(lam) {
    sum(
      evals * u_vec^2 /
        (evals + lam)^2
    )
  }
  
  norm_sq <- calc_norm_sq(lambda_min)
  
  if (norm_sq <= B_prime^2) {
    
    lambda_star <- lambda_min
    was_constrained <- FALSE
    
  } else {
    
    norm_func <- function(lam) {
      calc_norm_sq(lam) - B_prime^2
    }
    
    lam_upper <- max(1, lambda_min * 2)
    
    while (norm_func(lam_upper) > 0) {
      lam_upper <- lam_upper * 10
    }
    
    lambda_star <- uniroot(
      norm_func,
      lower = lambda_min,
      upper = lam_upper
    )$root
    
    was_constrained <- TRUE
  }
  
  beta <- as.vector(
    evecs %*%
      (
        u_vec / (evals + lambda_star)
      )
  )
  
  f_hat <- as.vector(K_T %*% beta)
  
  list(
    f_noisy_T = f_hat,
    lambda_used = lambda_star,
    was_constrained = was_constrained,
    norm_int = sqrt(max(norm_sq, 0)),
    norm_val = sqrt(
      max(calc_norm_sq(lambda_star), 0)
    )
  )
}

run_PTF_regime <- function(
    x_mat, design, B_prime, n, mu,
    lambda_min = 1e-8
) {
  
  prep <- prepare_PTF_data(
    x_mat = x_mat,
    B_prime = B_prime,
    n = n,
    mu = mu
  )
  
  fit <- fit_PTF_constrained(
    x_tilde = prep$x_tilde,
    design = design,
    B_prime = B_prime,
    lambda_min = lambda_min
  )
  
  c(
    list(
      x_tilde = prep$x_tilde,
      num_clamped = prep$num_clamped,
      prop_clamped = prep$prop_clamped,
      sigma_j = prep$sigma_j
    ),
    fit
  )
}

######## Mechanism 4: PTF with unconstrained KRR ########

fit_PTF_unconstrained <- function(
    x_tilde, design, lambda_min = 1e-8
) {
  
  K_T <- design$K_T
  T_pts <- design$T_pts
  
  beta <- as.vector(
    solve(
      K_T + diag(lambda_min, T_pts),
      x_tilde
    )
  )
  
  f_hat <- as.vector(K_T %*% beta)
  
  norm_sq <- as.numeric(
    crossprod(beta, K_T %*% beta)
  )
  
  list(
    f_noisy_T = f_hat,
    unconstrained_norm = sqrt(max(norm_sq, 0))
  )
}

run_PTF_regime_unconstrained <- function(
    x_mat, B_prime, n, mu,
    design, lambda_min = 1e-8
) {
  
  prep <- prepare_PTF_data(
    x_mat = x_mat,
    B_prime = B_prime,
    n = n,
    mu = mu
  )
  
  fit <- fit_PTF_unconstrained(
    x_tilde = prep$x_tilde,
    design = design,
    lambda_min = lambda_min
  )
  
  c(
    list(
      x_tilde = prep$x_tilde,
      num_clamped = prep$num_clamped,
      prop_clamped = prep$prop_clamped,
      sigma_j = prep$sigma_j
    ),
    fit
  )
}


######## Part 3. Simulation Function  

#### Run A1--A4 on one common noisy dataset ####
run_dp_mechanisms <- function(
    x_mat, design, B_prime, n, mu,
    lambda_min = 1e-8
) {
  
  # A1: FTP radial clamping
  A1 <- run_FTP_regime_clamp(
    x_mat = x_mat,
    design = design,
    B_prime = B_prime,
    n = n,
    mu = mu,
    lambda_min = lambda_min
  )
  
  # A3: FTP exact constrained fitting
  A3 <- run_FTP_regime_constrained(
    x_mat = x_mat,
    design = design,
    B_prime = B_prime,
    n = n,
    mu = mu,
    lambda_min = lambda_min
  )
  
  # One shared PTF preprocessing and DP realization
  PTF <- prepare_PTF_data(
    x_mat = x_mat,
    B_prime = B_prime,
    n = n,
    mu = mu
  )
  
  # A2 and A4 use the SAME x_tilde
  A2 <- fit_PTF_constrained(
    x_tilde = PTF$x_tilde,
    design = design,
    B_prime = B_prime,
    lambda_min = lambda_min
  )
  
  A4 <- fit_PTF_unconstrained(
    x_tilde = PTF$x_tilde,
    design = design,
    lambda_min = lambda_min
  )
  
  list(
    preds_FTP_clamp = A1$f_noisy_T,
    preds_PTF = A2$f_noisy_T,
    preds_FTP_constrained = A3$f_noisy_T,
    preds_PTF_unconstrained = A4$f_noisy_T,
    
    prop_clamp_FTP = A1$prop_clamped,
    prop_clamp_PTF = PTF$prop_clamped,
    
    num_clamp_FTP = A1$num_clamped,
    num_clamp_PTF = PTF$num_clamped,
    
    prop_constrained_FTP = A3$prop_constrained,
    num_constrained_FTP = A3$num_constrained,
    
    was_constrained_PTF = A2$was_constrained,
    norm_PTF = A2$norm_val,
    norm_PTF_int = A2$norm_int,
    lambda_star_PTF = A2$lambda_used,
    
    norms_FTP_clamp = A1$norms,
    norms_FTP_constrained = A3$unconstrained_norms,
    lambda_stars_FTP_constrained = A3$lambda_stars,
    
    norm_PTF_unconstrained = A4$unconstrained_norm
  )
}

#### Compute squared error ####
compute_ISE_RKHS <- function(
    f_bar, f_hat
) {
  
  if (length(f_bar) != length(f_hat)) {
    stop("f_bar and f_hat must have equal lengths.")
  }
  
  mean((f_hat - f_bar)^2)
}

######## One fixed RKHS simulation setting ########
run_one_setting_rkhs <- function(
    f_mat,
    f_bar,
    design,
    n,
    T_pts,
    M,
    B,
    B_prime,
    Bprime_choice,
    mu,
    sigma_0,
    sigma_kernel,
    lambda_min,
    reps
) {
  
  method_names <- c(
    "A1_FTP_clamp",
    "A2_PTF_constrained",
    "A3_FTP_constrained",
    "A4_PTF_unconstrained"
  )
  
  # Each row is one Monte Carlo repetition
  mse_mat <- matrix(
    NA_real_,
    nrow = reps,
    ncol = length(method_names),
    dimnames = list(NULL, method_names)
  )
  
  preds <- lapply(
    method_names,
    function(x) {
      matrix(NA_real_, nrow = reps, ncol = T_pts)
    }
  )
  
  names(preds) <- method_names
  
  # Diagnostics
  prop_clamp_FTP <- numeric(reps)
  prop_clamp_PTF <- numeric(reps)
  
  prop_constrained_FTP <- numeric(reps)
  prop_constrained_PTF <- numeric(reps)
  
  norm_PTF <- numeric(reps)
  norm_PTF_unconstrained <- numeric(reps)
  lambda_PTF <- numeric(reps)
  
  norm_FTP_clamp <- numeric(reps)
  norm_FTP_constrained <- numeric(reps)
  lambda_FTP_constrained <- numeric(reps)
  
  for (r in seq_len(reps)) {
    
    # 1. Fresh observation noise
    noisy <- generate_noisy_dataset(
      f_mat = f_mat,
      sigma_0 = sigma_0
    )
    
    # 2. All mechanisms share this noisy dataset
    out <- run_dp_mechanisms(
      x_mat = noisy$x_mat,
      design = design,
      B_prime = B_prime,
      n = n,
      mu = mu,
      lambda_min = lambda_min
    )
    
    current_preds <- list(
      A1_FTP_clamp = out$preds_FTP_clamp,
      A2_PTF_constrained = out$preds_PTF,
      A3_FTP_constrained = out$preds_FTP_constrained,
      A4_PTF_unconstrained = out$preds_PTF_unconstrained
    )
    
    # 3. Store predictions and repetition-level errors
    for (method in method_names) {
      
      preds[[method]][r, ] <- current_preds[[method]]
      
      mse_mat[r, method] <- compute_ISE_RKHS(
        f_bar = f_bar,
        f_hat = current_preds[[method]]
      )
    }
    
    # 4. Diagnostics
    prop_clamp_FTP[r] <- out$prop_clamp_FTP
    prop_clamp_PTF[r] <- out$prop_clamp_PTF
    
    prop_constrained_FTP[r] <-
      out$prop_constrained_FTP
    
    prop_constrained_PTF[r] <-
      as.numeric(out$was_constrained_PTF)
    
    norm_PTF[r] <- out$norm_PTF
    norm_PTF_unconstrained[r] <-
      out$norm_PTF_unconstrained
    
    lambda_PTF[r] <- out$lambda_star_PTF
    
    norm_FTP_clamp[r] <- mean(out$norms_FTP_clamp)
    
    norm_FTP_constrained[r] <-
      mean(out$norms_FTP_constrained)
    
    lambda_FTP_constrained[r] <-
      mean(out$lambda_stars_FTP_constrained)
  }
  
  ##########################################################
  ## Conditional IMSE and Monte Carlo SD
  ##########################################################
  
  IMSE <- colMeans(mse_mat)
  IMSE_sd <- apply(mse_mat, 2, sd)
  
  ##########################################################
  ## Empirical bias squared and variance
  ##
  ## All predictions target the fixed f_bar.
  ##########################################################
  
  bias2 <- numeric(length(method_names))
  variance <- numeric(length(method_names))
  
  names(bias2) <- method_names
  names(variance) <- method_names
  
  for (method in method_names) {
    
    pred_mat <- preds[[method]]
    
    pred_mean <- colMeans(pred_mat)
    
    bias2[method] <- mean(
      (pred_mean - f_bar)^2
    )
    
    # Population-normalized empirical variance
    pred_var <- apply(
      pred_mat, 2, var
    )
    
    variance[method] <- mean(
      pred_var * (reps - 1) / reps
    )
  }
  
  ##########################################################
  ## Theoretical DP-only variance: FTP
  ##########################################################
  
  sigma_DP <- 2 * B_prime / (n * mu)
  
  IVar_FTP <- sigma_DP^2 *
    mean(diag(design$K_T))
  
  ##########################################################
  ## Theoretical DP-only variance: PTF
  ##
  ## Exact for a fixed, non-data-dependent regularization
  ## parameter. For constrained A2, lambda is data-dependent,
  ## so this is not its exact unconditional DP variance.
  ##########################################################
  
  sigma_PTF <- 2 * B_prime * sqrt(T_pts) /
    (n * mu)
  
  S <- design$K_T %*% solve(
    design$K_T + diag(lambda_min, T_pts)
  )
  
  IVar_PTF_lambda_min <- sigma_PTF^2 *
    sum(S^2) / T_pts
  
  #### Ratios ####
  
  safe_ratio <- function(a, b) {
    if (is.finite(b) && b != 0) a / b else NA_real_
  }
  
  data.frame(
    n = n,
    T_pts = T_pts,
    M = M,
    B = B,
    mu = mu,
    sigma_0 = sigma_0,
    sigma_kernel = sigma_kernel,
    lambda_min = lambda_min,
    Bprime_choice = Bprime_choice,
    B_prime = B_prime,
    
    IMSE_A1_FTP_clamp = IMSE["A1_FTP_clamp"],
    IMSE_A2_PTF_constrained = IMSE["A2_PTF_constrained"],
    IMSE_A3_FTP_constrained = IMSE["A3_FTP_constrained"],
    IMSE_A4_PTF_unconstrained = IMSE["A4_PTF_unconstrained"],
    
    IMSE_A1_sd = IMSE_sd["A1_FTP_clamp"],
    IMSE_A2_sd = IMSE_sd["A2_PTF_constrained"],
    IMSE_A3_sd = IMSE_sd["A3_FTP_constrained"],
    IMSE_A4_sd = IMSE_sd["A4_PTF_unconstrained"],
    
    Bias2_A1 = bias2["A1_FTP_clamp"],
    Bias2_A2 = bias2["A2_PTF_constrained"],
    Bias2_A3 = bias2["A3_FTP_constrained"],
    Bias2_A4 = bias2["A4_PTF_unconstrained"],
    
    Var_A1 = variance["A1_FTP_clamp"],
    Var_A2 = variance["A2_PTF_constrained"],
    Var_A3 = variance["A3_FTP_constrained"],
    Var_A4 = variance["A4_PTF_unconstrained"],
    
    IVar_FTP = IVar_FTP,
    IVar_PTF_lambda_min = IVar_PTF_lambda_min,
    
    Ratio_A1_A2 = safe_ratio(IMSE[1], IMSE[2]),
    Ratio_A3_A2 = safe_ratio(IMSE[3], IMSE[2]),
    Ratio_A1_A3 = safe_ratio(IMSE[1], IMSE[3]),
    Ratio_A2_A4 = safe_ratio(IMSE[2], IMSE[4]),
    
    A1_prop_clamped = mean(prop_clamp_FTP),
    PTF_prop_clamped = mean(prop_clamp_PTF),
    A3_prop_constrained = mean(prop_constrained_FTP),
    A2_prop_constrained = mean(prop_constrained_PTF),
    
    FTP_norm_clamp_mean = mean(norm_FTP_clamp),
    FTP_norm_constrained_mean = mean(norm_FTP_constrained),
    FTP_lambda_constrained_mean = mean(lambda_FTP_constrained),
    
    PTF_norm_mean = mean(norm_PTF),
    PTF_lambda_mean = mean(lambda_PTF),
    PTF_unconstrained_norm_mean = mean(norm_PTF_unconstrained),
    
    stringsAsFactors = FALSE
  )
}

######## Main RKHS simulation ########
run_simulation_rkhs <- function(
    
  n_vec = c(50, 100, 200, 400),
  
  mu_vec = c(0.1, 0.5, 1),
  
  T_pts_vec = c(5, 10, 20),
  
  M_vec = c(10, 20),
  
  B_vec = c(0.5, 1, 2),
  
  Bprime_choices = c(
    "B", "sqrtlog", "log", "sqrtn", "n_2_3"
  ),
  
  sigma_0_vec = c(0, 0.1, 1),
  
  sigma_kernel_vec = c(0.2),
  
  lambda_min_vec = c(1e-8),
  
  mu0 = NULL,
  
  reps = 1000,
  
  n_cores = NULL,
  
  seed = 12345
) {
  
  ##########################################################
  ## 1. Validate parameters
  ##########################################################
  
  if (any(n_vec < 1) ||
      any(T_pts_vec < 1) ||
      any(M_vec < 1)) {
    stop("n_vec, T_pts_vec, and M_vec must be positive.")
  }
  
  if (any(mu_vec <= 0)) {
    stop("All privacy parameters mu must be positive.")
  }
  
  if (any(B_vec <= 0)) {
    stop("All B values must be positive.")
  }
  
  if (any(sigma_kernel_vec <= 0)) {
    stop("All kernel bandwidths must be positive.")
  }
  
  if (any(lambda_min_vec <= 0)) {
    stop("lambda_min_vec must contain strictly positive values.")
  }
  
  if (any(sigma_0_vec < 0)) {
    stop("Observation noise standard deviations must be nonnegative.")
  }
  
  if (reps < 2) {
    stop("Use at least two repetitions to estimate variance.")
  }
  
  ##########################################################
  ## 2. Generate the complete parameter grid
  ##########################################################
  
  grid <- expand.grid(
    n = n_vec,
    T_pts = T_pts_vec,
    M = M_vec,
    B = B_vec,
    mu = mu_vec,
    sigma_0 = sigma_0_vec,
    sigma_kernel = sigma_kernel_vec,
    lambda_min = lambda_min_vec,
    Bprime_choice = Bprime_choices,
    stringsAsFactors = FALSE
  )
  
  cat(
    "Total parameter settings:",
    nrow(grid), "\n"
  )
  
  ##########################################################
  ## 3. Fixed latent datasets
  ##########################################################
  
  set.seed(seed)
  
  latent_data <- list()
  
  latent_grid <- unique(
    grid[
      ,
      c("n", "M", "B", "sigma_kernel"),
      drop = FALSE
    ]
  )
  
  for (i in seq_len(nrow(latent_grid))) {
    
    n_current <- latent_grid$n[i]
    M_current <- latent_grid$M[i]
    B_current <- latent_grid$B[i]
    sk_current <- latent_grid$sigma_kernel[i]
    
    key <- paste(
      n_current,
      M_current,
      B_current,
      sk_current,
      sep = "_"
    )
    
    # Default mu0 varies with M, as in the original RKHS code.
    # If explicitly supplied, it must have length M.
    current_mu0 <- mu0
    
    if (!is.null(current_mu0) &&
        length(current_mu0) != M_current) {
      stop(
        "When M varies, mu0 must be NULL or have length M ",
        "for every selected M."
      )
    }
    
    cat(
      "Generating fixed latent data: n =", n_current,
      ", M =", M_current,
      ", B =", B_current,
      ", sigma_kernel =", sk_current, "\n"
    )
    
    latent_data[[key]] <- generate_rkhs_latent(
      n = n_current,
      M = M_current,
      B = B_current,
      sigma_kernel = sk_current,
      mu0 = current_mu0
    )
  }
  
  ##########################################################
  ## 4. Prepare fixed kernel designs
  ##########################################################
  
  design_data <- list()
  
  design_grid <- unique(
    grid[
      ,
      c("T_pts", "sigma_kernel"),
      drop = FALSE
    ]
  )
  
  for (i in seq_len(nrow(design_grid))) {
    
    T_current <- design_grid$T_pts[i]
    sk_current <- design_grid$sigma_kernel[i]
    
    key <- paste(
      T_current,
      sk_current,
      sep = "_"
    )
    
    cat(
      "Preparing kernel design: T_pts =", T_current,
      ", sigma_kernel =", sk_current, "\n"
    )
    
    design_data[[key]] <- prepare_rkhs_design(
      T_pts = T_current,
      sigma_kernel = sk_current
    )
  }
  
  ##########################################################
  ## 5. Configure parallel workers
  ##########################################################
  
  if (is.null(n_cores)) {
    n_cores <- max(
      1,
      parallel::detectCores() - 2
    )
  }
  
  n_cores <- min(
    n_cores,
    nrow(grid),
    100
  )
  
  cat(
    "Using", n_cores,
    "parallel workers for",
    nrow(grid), "settings.\n"
  )
  
  cl <- parallel::makeCluster(
    n_cores,
    type = "PSOCK"
  )
  
  on.exit(
    parallel::stopCluster(cl),
    add = TRUE
  )
  
  parallel::clusterSetRNGStream(
    cl,
    iseed = seed + 1
  )
  
  doParallel::registerDoParallel(cl)
  
  ##########################################################
  ## 6. Run parameter settings in parallel
  ##########################################################
  
  results <- foreach(
    k = seq_len(nrow(grid)),
    .combine = rbind,
    .packages = character(0),
    .export = c(
      "rbf_kernel",
      "get_Bprime",
      "generate_noisy_dataset",
      "prepare_rkhs_design",
      "evaluate_rkhs_data",
      "run_FTP_regime_clamp",
      "run_FTP_regime_constrained",
      "prepare_PTF_data",
      "fit_PTF_constrained",
      "fit_PTF_unconstrained",
      "run_dp_mechanisms",
      "compute_ISE_RKHS",
      "run_one_setting_rkhs"
    )
  ) %dopar% {
    
    ########################################################
    ## Current setting
    ########################################################
    
    n_current <- grid$n[k]
    T_current <- grid$T_pts[k]
    M_current <- grid$M[k]
    B_current <- grid$B[k]
    mu_current <- grid$mu[k]
    sigma_current <- grid$sigma_0[k]
    sk_current <- grid$sigma_kernel[k]
    lambda_current <- grid$lambda_min[k]
    Bprime_choice <- grid$Bprime_choice[k]
    
    ########################################################
    ## Retrieve fixed latent data
    ########################################################
    
    latent_key <- paste(
      n_current,
      M_current,
      B_current,
      sk_current,
      sep = "_"
    )
    
    latent <- latent_data[[latent_key]]
    
    # Evaluate the same latent coefficient vectors on this grid
    evaluated <- evaluate_rkhs_data(
      alpha_mat = latent$alpha_mat,
      T_pts = T_current,
      sigma_kernel = sk_current
    )
    
    ########################################################
    ## Retrieve fixed kernel design
    ########################################################
    
    design_key <- paste(
      T_current,
      sk_current,
      sep = "_"
    )
    
    design <- design_data[[design_key]]
    
    ########################################################
    ## B-prime
    ########################################################
    
    B_prime <- get_Bprime(
      n = n_current,
      B = B_current,
      choice = Bprime_choice
    )
    
    ########################################################
    ## One complete Monte Carlo setting
    ########################################################
    
    run_one_setting_rkhs(
      f_mat = evaluated$f_mat,
      f_bar = evaluated$f_bar,
      design = design,
      
      n = n_current,
      T_pts = T_current,
      M = M_current,
      
      B = B_current,
      B_prime = B_prime,
      Bprime_choice = Bprime_choice,
      
      mu = mu_current,
      sigma_0 = sigma_current,
      sigma_kernel = sk_current,
      lambda_min = lambda_current,
      
      reps = reps
    )
  }
  
  ##########################################################
  ## 7. Sort and return
  ##########################################################
  
  rownames(results) <- NULL
  
  results <- results[
    order(
      results$n,
      results$M,
      results$T_pts,
      results$B,
      results$sigma_0,
      results$sigma_kernel,
      results$lambda_min,
      results$mu,
      match(
        results$Bprime_choice,
        c("B", "sqrtlog", "log", "sqrtn", "n_2_3")
      )
    ),
  ]
  
  rownames(results) <- NULL
  
  cat("RKHS simulation complete.\n")
  
  results
}

##### Simulation 1 for Noiseless

#tbl_all <- run_simulation(
# n_vec = 2^(1:9),
# mu_vec = c(0.1, 1),
# T_pts_vec = c(2, 3, 4, 5, 6, 8, 10),
# M_vec = c(20),
# B_vec = c(1),
# Bprime_choices = c("sqrtlog"),
# sigma_0_vec = c(0),
# sigma_kernel_vec = c(0.2),
# lambda_min_vec = c(1e-8),
# mu0 = NULL,
# reps = 5000,
# n_cores = NULL,
# seed = 12345
#)

##### Simulation 2 for B' choices
# tbl_all <- run_simulation(
# n_vec = 2^(2:11),
# mu_vec = c(1),
# T_pts_vec = c(4),
# M_vec = c(20),
# B_vec = c(1),
# Bprime_choices = c("B", "sqrtlog", "log", "sqrtn", "n_2_3"),
# sigma_0_vec = c(1),
# sigma_kernel_vec = c(0.2),
# lambda_min_vec = c(1e-8),
# mu0 = NULL,
# reps = 5000,
# n_cores = NULL,
# seed = 12345
# )

##### Simulation 3 for Noisy Level
#tbl_all <- run_simulation(
# n_vec = 2^(2:11),
# mu_vec = c(1),
# T_pts_vec = c(2, 3, 4, 5, 6, 8, 10),
# M_vec = c(20),
# B_vec = c(1),
# Bprime_choices = c("sqrtlog"),
# sigma_0_vec = c(0.01, 1, 100),
# sigma_kernel_vec = c(0.2),
# lambda_min_vec = c(1e-8),
# mu0 = NULL,
# reps = 5000,
# n_cores = NULL,
# seed = 12345
#)

#### Run RKHS simulation ####
tbl_all <- run_simulation_rkhs(
  n_vec = 2^(2:11),
  mu_vec = c(1),
  T_pts_vec = c(2, 3, 4, 5, 6, 8, 10),
  M_vec = c(20),
  B_vec = c(1),
  Bprime_choices = c("log"),
  sigma_0_vec = c(0),
  sigma_kernel_vec = c(0.2),
  lambda_min_vec = c(1e-8),
  mu0 = NULL,
  reps = 5000,
  n_cores = NULL,
  seed = 12345
)

tbl_all <- tbl_all[
  order(
    tbl_all$n,
    tbl_all$M,
    tbl_all$T_pts,
    tbl_all$B,
    tbl_all$sigma_0,
    tbl_all$sigma_kernel,
    tbl_all$lambda_min,
    tbl_all$mu,
    match(
      tbl_all$Bprime_choice,
      c("B", "sqrtlog", "log", "sqrtn", "n_2_3")
    )
  ),
]

rownames(tbl_all) <- NULL

# Save it
saveRDS(tbl_all, file = "RKHS_sigma0_0_mu_1_n_4-2048.rds")
