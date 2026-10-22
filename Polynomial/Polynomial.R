############################################################
## Polynomial Case:
## Fixed theta / Fixed latent functions,
## Random observation noise + Random DP noise
##
## A1: FTP + radial clamping
## A2: PTF + constrained WLS
## A3: FTP + exact constrained least squares
## A4: PTF + unconstrained WLS
##
## For each n:
##   1. Generate ONE theta dataset.
##   2. Keep theta fixed for all MC repetitions.
##   3. Compute noiseless functions f_i = Phi theta_i once.
##   4. For every MC repetition:
##        - generate fresh observation noise
##        - construct noisy observations
##        - perform deterministic preprocessing
##        - generate fresh DP noise
##        - estimate A1--A4
##
## Thus the Monte Carlo randomness comes from:
##   (i) observation noise
##   (ii) DP noise
##
## Quantities depending only on the fixed design are computed once.
##
## B' choices:
##   B
##   sqrt(log n)
##   log(n)
##   sqrt(n)
##   n^(2/3)
############################################################

library("doParallel")
library("foreach")

############################################################
##### General shifted orthonormal Legendre basis
############################################################

legendre_poly <- function(k, x){
  if(k==0) return(rep(1,length(x)))
  if(k==1) return(x)
  Pkm1 <- rep(1,length(x))
  Pk <- x
  for(i in 1:(k-1)){
    Pkp1 <- ((2*i+1)*x*Pk - i*Pkm1)/(i+1)
    Pkm1 <- Pk
    Pk <- Pkp1
  }
  Pk
}

phi_basis <- function(t, m) {
  x <- 2 * t - 1
  Phi <- matrix(0, nrow = length(t), ncol = m)
  for (k in seq_len(m)) {
    Phi[, k] <- sqrt(2 * k - 1) * legendre_poly(k - 1, x)
  }
  Phi
}

############################################################
##### Fixed latent coefficient generation
############################################################

generate_fixed_theta <- function(n, m, B, mu0 = rep(0, m)) {
  
  # Gaussian directions: Z_i ~ N(mu, I_m)
  Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
  Z <- sweep(Z, 2, mu0, "+")
  
  # Normalize each observation to unit l2 norm
  Z_norm <- sqrt(rowSums(Z^2))
  U <- Z / Z_norm
  
  # Uniform in volume of m-dimensional ball
  r <- B * runif(n)^(1 / m)
  
  # Scale each row
  theta <- U * r
  
  list(
    theta = theta,
    theta_bar_true = colMeans(theta)
  )
}

############################################################
##### Function evaluation
############################################################

compute_f_matrix <- function(theta, Phi) {
  theta %*% t(Phi)
}

############################################################
##### B' choices
############################################################

get_Bprime <- function(n, B, choice) {
  
  switch(
    choice,
    B = B,
    sqrtlog = sqrt(log(n)),
    log = log(n),
    sqrtn =sqrt(n),
    n_2_3 = n^(2/3),
    stop("Unknown Bprime choice.")
  )
}

############################################################
##### Constrained quadratic solver
############################################################

solve_constrained_from_eigen <- function(
    c_vec,
    d,
    U,
    Bprime,
    tol = 1e-10
) {
  
  ## Unconstrained solution in eigenbasis:
  ##
  ## theta_OLS = U diag(1/d) U' c
  
  theta_ols <- as.vector(
    U %*% (c_vec / d)
  )
  
  norm_ols <- sqrt(
    sum(theta_ols^2)
  )
  
  ## Constraint inactive
  if (norm_ols <= Bprime) {
    
    return(
      list(
        theta = theta_ols,
        was_constrained = FALSE
      )
    )
  }
  
  ## Constraint active.
  ##
  ## Solve:
  ##
  ## sum_j c_j^2 / (d_j + lambda)^2 = Bprime^2
  
  norm_equation <- function(lambda) {
    
    sum(
      c_vec^2 /
        (d + lambda)^2
    ) - Bprime^2
  }
  
  lambda_upper <- 1
  
  while (
    norm_equation(lambda_upper) > 0
  ) {
    lambda_upper <- 2 * lambda_upper
  }
  
  lambda_star <- uniroot(
    norm_equation,
    lower = 0,
    upper = lambda_upper,
    tol = tol
  )$root
  
  theta_star <- as.vector(
    U %*%
      (
        c_vec /
          (d + lambda_star)
      )
  )
  
  list(
    theta = theta_star,
    was_constrained = TRUE
  )
}


############################################################
##### Fixed-design linear algebra
############################################################

prepare_design <- function(Phi) {
  
  XtX <- crossprod(Phi)
  
  XtX <- (
    XtX + t(XtX)
  ) / 2
  
  eig <- eigen(
    XtX,
    symmetric = TRUE
  )
  
  if (
    any(eig$values <= 1e-12)
  ) {
    stop(
      "Polynomial design matrix is numerically rank deficient."
    )
  }
  
  list(
    XtX = XtX,
    XtX_inv = solve(XtX),
    d = eig$values,
    U = eig$vectors
  )
}


############################################################
##### Generate noisy observations
############################################################
##
## theta and f_mat are fixed.
## Only observation noise changes across MC repetitions.
############################################################

generate_noisy_dataset <- function(
    f_mat,
    sigma_obs,
    Phi,
    design
) {
  
  n <- nrow(f_mat)
  T <- ncol(f_mat)
  
  ##########################################################
  ## Fresh observation noise
  ##########################################################
  
  f_mat_noisy <- f_mat +
    matrix(
      rnorm(
        n * T,
        mean = 0,
        sd = sigma_obs
      ),
      nrow = n,
      ncol = T
    )
  
  ##########################################################
  ## Individual unconstrained OLS coefficients
  ##########################################################
  
  XtY <- crossprod(
    Phi,
    t(f_mat_noisy)
  )
  
  theta_ols <- t(
    solve(
      design$XtX,
      XtY
    )
  )
  
  list(
    f_mat_noisy = f_mat_noisy,
    theta_ols = theta_ols
  )
}


############################################################
##### A1: deterministic preprocessing
############################################################

prepare_A1 <- function(
    theta_ols,
    Bprime
) {
  
  norms <- sqrt(
    rowSums(theta_ols^2)
  )
  
  theta_protected <- theta_ols
  
  idx <- norms > Bprime
  
  if (any(idx)) {
    
    theta_protected[idx, ] <-
      theta_ols[idx, , drop = FALSE] *
      (
        Bprime /
          norms[idx]
      )
  }
  
  list(
    theta_mean = colMeans(
      theta_protected
    ),
    prop_clamped = mean(idx)
  )
}


############################################################
##### A1: DP randomization
############################################################

privatize_A1 <- function(
    prep,
    n,
    Bprime,
    mu
) {
  
  Delta <- 2 * Bprime / n
  
  sigma_DP <- Delta / mu
  
  list(
    theta =
      prep$theta_mean +
      rnorm(
        length(prep$theta_mean),
        mean = 0,
        sd = sigma_DP
      ),
    
    Delta = Delta,
    sigma_DP = sigma_DP
  )
}


############################################################
##### A3: deterministic preprocessing
############################################################

prepare_A3 <- function(
    f_mat_noisy,
    Phi,
    design,
    Bprime
) {
  
  ##########################################################
  ## Phi' y_i for every individual
  ##########################################################
  
  XtY <- crossprod(
    Phi,
    t(f_mat_noisy)
  )
  
  ##########################################################
  ## Transform linear terms into eigenbasis
  ##########################################################
  
  C <- crossprod(
    design$U,
    XtY
  )
  
  n <- nrow(
    f_mat_noisy
  )
  
  theta_hat <- matrix(
    0,
    nrow = n,
    ncol = ncol(Phi)
  )
  
  was_constrained <- logical(n)
  
  ##########################################################
  ## Solve constrained problem individually
  ##########################################################
  
  for (i in seq_len(n)) {
    
    sol <- solve_constrained_from_eigen(
      c_vec = C[, i],
      d = design$d,
      U = design$U,
      Bprime = Bprime
    )
    
    theta_hat[i, ] <- sol$theta
    
    was_constrained[i] <-
      sol$was_constrained
  }
  
  list(
    theta_mean = colMeans(
      theta_hat
    ),
    
    prop_constrained =
      mean(was_constrained)
  )
}


############################################################
##### A3: DP randomization
############################################################

privatize_A3 <- function(
    prep,
    n,
    Bprime,
    mu
) {
  
  Delta <- 2 * Bprime / n
  
  sigma_DP <- Delta / mu
  
  list(
    theta =
      prep$theta_mean +
      rnorm(
        length(prep$theta_mean),
        mean = 0,
        sd = sigma_DP
      ),
    
    Delta = Delta,
    sigma_DP = sigma_DP
  )
}


############################################################
##### PTF: deterministic preprocessing
############################################################

prepare_PTF <- function(
    f_mat_noisy,
    Phi,
    n,
    Bprime,
    phi_norm
) {
  
  ##########################################################
  ## Individual coordinate bounds:
  ##
  ## |f_i(t_j)| <= B' ||phi(t_j)||
  ##########################################################
  
  upper <- Bprime * phi_norm
  
  upper_mat <- matrix(
    upper,
    nrow = n,
    ncol = ncol(Phi),
    byrow = TRUE
  )
  
  ##########################################################
  ## Coordinate-wise clipping
  ##########################################################
  
  f_clip <- pmin(
    pmax(
      f_mat_noisy,
      -upper_mat
    ),
    upper_mat
  )
  
  clipped <- (
    abs(
      f_mat_noisy - f_clip
    ) > 0
  )
  
  ##########################################################
  ## Average after clipping
  ##########################################################
  
  xbar <- colMeans(
    f_clip
  )
  
  ##########################################################
  ## Sensitivity
  ##########################################################
  
  Delta_vec <-
    2 * Bprime * phi_norm / n
  
  list(
    xbar = xbar,
    Delta_vec = Delta_vec,
    prop_clamped = mean(clipped)
  )
}


############################################################
##### PTF: DP randomization
############################################################

privatize_PTF <- function(
    prep,
    T,
    mu
) {
  
  sigma_vec <-
    prep$Delta_vec *
    sqrt(T) /
    mu
  
  xbar_noisy <-
    prep$xbar +
    rnorm(
      T,
      mean = 0,
      sd = sigma_vec
    )
  
  list(
    xbar_noisy = xbar_noisy,
    sigma_vec = sigma_vec
  )
}


############################################################
##### PTF WLS design
############################################################
##
## The WLS weights have the form
##
##   1 / sigma_j^2
##
## and the common scalar does not affect the WLS point
## estimate.
############################################################

prepare_PTF_WLS <- function(
    Phi,
    phi_norm
) {
  
  w_shape <- 1 / phi_norm^2
  
  XtWX <- crossprod(
    Phi,
    w_shape * Phi
  )
  
  XtWX <- (
    XtWX + t(XtWX)
  ) / 2
  
  eig <- eigen(
    XtWX,
    symmetric = TRUE
  )
  
  if (
    any(eig$values <= 1e-12)
  ) {
    stop(
      "PTF WLS design matrix is numerically rank deficient."
    )
  }
  
  list(
    weights_shape = w_shape,
    XtWX = XtWX,
    XtWX_inv = solve(XtWX),
    d = eig$values,
    U = eig$vectors
  )
}


############################################################
##### A2: constrained WLS
############################################################

fit_A2 <- function(
    xbar_noisy,
    Phi,
    PTF_WLS,
    Bprime
) {
  
  XtWy <- crossprod(
    Phi,
    PTF_WLS$weights_shape *
      xbar_noisy
  )
  
  c_vec <- crossprod(
    PTF_WLS$U,
    XtWy
  )
  
  solve_constrained_from_eigen(
    c_vec = c_vec,
    d = PTF_WLS$d,
    U = PTF_WLS$U,
    Bprime = Bprime
  )
}


############################################################
##### A4: unconstrained WLS
############################################################

fit_A4 <- function(
    xbar_noisy,
    Phi,
    PTF_WLS
) {
  
  XtWy <- crossprod(
    Phi,
    PTF_WLS$weights_shape *
      xbar_noisy
  )
  
  theta <- as.vector(
    solve(
      PTF_WLS$XtWX,
      XtWy
    )
  )
  
  list(
    theta = theta
  )
}


############################################################
##### ISE
############################################################

compute_ISE <- function(
    theta_bar_true,
    theta_hat
) {
  
  sum(
    (
      theta_hat -
        theta_bar_true
    )^2
  )
}


############################################################
##### Theoretical DP variance: FTP
############################################################

compute_IVar_FTP <- function(
    Bprime,
    n,
    mu,
    m
) {
  
  sigma_DP <-
    2 * Bprime /
    (n * mu)
  
  m *
    sigma_DP^2
}


############################################################
##### Theoretical DP variance: PTF
############################################################

compute_IVar_PTF <- function(
    Bprime,
    n,
    mu,
    T,
    Phi,
    phi_norm
) {
  
  w_shape <-
    1 / phi_norm^2
  
  XtWX_shape <-
    crossprod(
      Phi,
      w_shape * Phi
    )
  
  ##########################################################
  ## Full WLS weights:
  ##
  ## w_j =
  ## n^2 mu^2 /
  ## (4 B'^2 T)
  ## * 1 / ||phi(t_j)||^2
  ##
  ## Hence
  ##
  ## Cov(theta)
  ## =
  ## 4 B'^2 T /
  ## (n^2 mu^2)
  ## * (XtWX_shape)^(-1)
  ##########################################################
  
  scale_cov <-
    4 * Bprime^2 * T /
    (n^2 * mu^2)
  
  Cov_theta <-
    scale_cov *
    solve(XtWX_shape)
  
  sum(
    diag(Cov_theta)
  )
}

############################################################
##### One (n, B', mu) simulation
############################################################

run_one_setting <- function(
    theta,
    f_mat,
    theta_bar_true,
    Phi,
    design,
    PTF_WLS,
    n,
    T,
    m,
    mu,
    B,
    Bprime,
    Bprime_choice,
    sigma_obs,
    phi_norm,
    reps
) {
  
  ##########################################################
  ## Storage
  ##########################################################
  
  ISE_A1 <- numeric(reps)
  ISE_A2 <- numeric(reps)
  ISE_A3 <- numeric(reps)
  ISE_A4 <- numeric(reps)
  
  A1_clamped_prop <- numeric(reps)
  A3_constrained_prop <- numeric(reps)
  PTF_clipped_prop <- numeric(reps)
  A2_constrained <- logical(reps)
  
  ##########################################################
  ## Monte Carlo loop
  ##########################################################
  
  for (r in seq_len(reps)) {
    
    ########################################################
    ## 1. Fresh observation noise
    ########################################################
    
    noisy <- generate_noisy_dataset(
      f_mat = f_mat,
      sigma_obs = sigma_obs,
      Phi = Phi,
      design = design
    )
    
    ########################################################
    ## 2. Deterministic preprocessing conditional on
    ##    this observation-noise realization
    ########################################################
    
    A1 <- prepare_A1(
      theta_ols = noisy$theta_ols,
      Bprime = Bprime
    )
    
    A3 <- prepare_A3(
      f_mat_noisy = noisy$f_mat_noisy,
      Phi = Phi,
      design = design,
      Bprime = Bprime
    )
    
    PTF <- prepare_PTF(
      f_mat_noisy = noisy$f_mat_noisy,
      Phi = Phi,
      n = n,
      Bprime = Bprime,
      phi_norm = phi_norm
    )
    
    ########################################################
    ## 3. Fresh DP randomness
    ########################################################
    
    A1_DP <- privatize_A1(
      prep = A1,
      n = n,
      Bprime = Bprime,
      mu = mu
    )
    
    A3_DP <- privatize_A3(
      prep = A3,
      n = n,
      Bprime = Bprime,
      mu = mu
    )
    
    ########################################################
    ## Same privatized PTF data for A2 and A4
    ########################################################
    
    PTF_DP <- privatize_PTF(
      prep = PTF,
      T = T,
      mu = mu
    )
    
    ########################################################
    ## 4. PTF estimators
    ########################################################
    
    A2_fit <- fit_A2(
      xbar_noisy = PTF_DP$xbar_noisy,
      Phi = Phi,
      PTF_WLS = PTF_WLS,
      Bprime = Bprime
    )
    
    A4_fit <- fit_A4(
      xbar_noisy = PTF_DP$xbar_noisy,
      Phi = Phi,
      PTF_WLS = PTF_WLS
    )
    
    ########################################################
    ## 5. Store ISE
    ########################################################
    
    ISE_A1[r] <- compute_ISE(
      theta_bar_true = theta_bar_true,
      theta_hat = A1_DP$theta
    )
    
    ISE_A2[r] <- compute_ISE(
      theta_bar_true = theta_bar_true,
      theta_hat = A2_fit$theta
    )
    
    ISE_A3[r] <- compute_ISE(
      theta_bar_true = theta_bar_true,
      theta_hat = A3_DP$theta
    )
    
    ISE_A4[r] <- compute_ISE(
      theta_bar_true = theta_bar_true,
      theta_hat = A4_fit$theta
    )
    
    ########################################################
    ## 6. Store preprocessing diagnostics
    ########################################################
    
    A1_clamped_prop[r] <-
      A1$prop_clamped
    
    A3_constrained_prop[r] <-
      A3$prop_constrained
    
    PTF_clipped_prop[r] <-
      PTF$prop_clamped
    
    A2_constrained[r] <-
      A2_fit$was_constrained
  }
  
  
  ##########################################################
  ## Theoretical DP-only variance
  ##
  ## These depend only on:
  ## B', n, mu, T, Phi
  ##
  ## They do not need MC averaging.
  ##########################################################
  
  IVar_FTP <- compute_IVar_FTP(
    Bprime = Bprime,
    n = n,
    mu = mu,
    m = m
  )
  
  IVar_PTF <- compute_IVar_PTF(
    Bprime = Bprime,
    n = n,
    mu = mu,
    T = T,
    Phi = Phi,
    phi_norm = phi_norm
  )
  
  
  ##########################################################
  ## Return results
  ##########################################################
  
  data.frame(
    
    n = n,
    T = T,
    m = m,
    mu = mu,
    B = B,
    sigma_obs = sigma_obs,
    Bprime_choice =
      Bprime_choice,
    
    Bprime =
      Bprime,
    
    
    ########################################################
    ## IMSE
    ########################################################
    
    IMSE_A1_FTP_clamp =
      mean(ISE_A1),
    
    IMSE_A2_PTF_constrained =
      mean(ISE_A2),
    
    IMSE_A3_FTP_constrained =
      mean(ISE_A3),
    
    IMSE_A4_PTF_unconstrained =
      mean(ISE_A4),
    
    
    ########################################################
    ## Monte Carlo SD
    ########################################################
    
    IMSE_A1_sd =
      sd(ISE_A1),
    
    IMSE_A2_sd =
      sd(ISE_A2),
    
    IMSE_A3_sd =
      sd(ISE_A3),
    
    IMSE_A4_sd =
      sd(ISE_A4),
    
    
    ########################################################
    ## DP-only variance
    ########################################################
    
    IVar_FTP =
      IVar_FTP,
    
    IVar_PTF =
      IVar_PTF,
    
    
    ########################################################
    ## Ratios
    ########################################################
    
    Ratio_A1_A2 =
      mean(ISE_A1) /
      mean(ISE_A2),
    
    Ratio_A3_A2 =
      mean(ISE_A3) /
      mean(ISE_A2),
    
    Ratio_A1_A3 =
      mean(ISE_A1) /
      mean(ISE_A3),
    
    Ratio_A2_A4 =
      mean(ISE_A2) /
      mean(ISE_A4),
    
    
    ########################################################
    ## Average preprocessing activity
    ########################################################
    
    A1_prop_clamped =
      mean(A1_clamped_prop),
    
    A2_prop_constrained =
      mean(A2_constrained),
    
    A3_prop_constrained =
      mean(A3_constrained_prop),
    
    PTF_prop_clamped =
      mean(PTF_clipped_prop)
  )
}

############################################################
##### Main simulation
############################################################

run_simulation <- function(
    n_vec = c(
      50,
      100,
      200,
      400,
      800,
      1600,
      3200
    ),
    
    mu_vec =
      seq(
        0.1,
        1.5,
        by = 0.2
      ),
    
    Bprime_choices =
      c(
        "B",
        "sqrtlog",
        "log",
        "sqrtn",
        "n_2_3"
      ),
    
    T_vec = c(20, 30, 50),
    m_vec = c(2, 3, 4, 5, 6, 8, 10),
    B_vec = c(0.5, 1, 2),
    mu0 = NULL,
    
    sigma_obs_vec = c(0, 0.1, 1, 10),
    
    reps = 1000,
    
    n_cores = NULL,
    
    seed = 12345
) {
  
  
  ############################################################
  ## Fixed designs for every (T, m)
  ############################################################
  
  if (any(T_vec < max(m_vec))) {
    stop("Every T must be at least max(m_vec).")
  }
  
  design_data <- list()
  
  for (T_current in T_vec) {
    
    t <- seq(0, 1, length.out = T_current)
    
    for (m_current in m_vec) {
      
      cat(
        "Preparing design: T =", T_current,
        ", m =", m_current, "\n"
      )
      
      Phi <- phi_basis(t, m_current)
      
      design <- prepare_design(Phi)
      
      phi_norm <- sqrt(rowSums(Phi^2))
      
      PTF_WLS <- prepare_PTF_WLS(
        Phi = Phi,
        phi_norm = phi_norm
      )
      
      key <- paste(T_current, m_current, sep = "_")
      
      design_data[[key]] <- list(
        Phi = Phi,
        design = design,
        phi_norm = phi_norm,
        PTF_WLS = PTF_WLS
      )
    }
  }
  
  ############################################################
  ## One fixed theta dataset per (n, m, B)
  ############################################################
  
  set.seed(seed)
  
  fixed_data <- list()
  
  for (m_current in m_vec) {
    
    for (B_current in B_vec) {
      
      mu0_current <- if (is.null(mu0)) {
        seq(-1, 1, length.out = m_current)
      } else {
        if (length(mu0) != m_current) {
          stop(
            "mu0 must have length m, or be NULL ",
            "when m varies."
          )
        }
        mu0
      }
      
      for (n in n_vec) {
        
        cat(
          "Generating fixed theta: n =", n,
          ", m =", m_current,
          ", B =", B_current, "\n"
        )
        
        truth <- generate_fixed_theta(
          n = n,
          m = m_current,
          B = B_current,
          mu0 = mu0_current
        )
        
        key <- paste(n, m_current, B_current, sep = "_")
        
        fixed_data[[key]] <- list(
          theta = truth$theta,
          theta_bar_true = truth$theta_bar_true
        )
      }
    }
  }
  
  
  ##########################################################
  ## Parallel task grid
  ##########################################################
  
  grid <- expand.grid(
    n = n_vec,
    m = m_vec,
    T = T_vec,
    B = B_vec,
    sigma_obs = sigma_obs_vec,
    Bprime_choice = Bprime_choices,
    mu = mu_vec,
    stringsAsFactors = FALSE
  )
  
  
  ##########################################################
  ## Number of workers
  ##########################################################
  
  if (is.null(n_cores)) {
    
    n_cores <- min(
      
      max(
        1,
        parallel::detectCores() - 2
      ),
      
      nrow(grid)
    )
  }
  
  n_cores <-
    min(
      n_cores,
      nrow(grid),
      100
    )
  
  cat(
    "Using",
    n_cores,
    "parallel workers for",
    nrow(grid),
    "settings.\n"
  )
  
  
  ##########################################################
  ## Parallel cluster
  ##########################################################
  
  cl <- parallel::makeCluster(
    n_cores,
    type = "PSOCK"
  )
  
  on.exit(
    parallel::stopCluster(cl),
    add = TRUE
  )
  
  ##########################################################
  ## Independent reproducible RNG streams
  ##########################################################
  
  parallel::clusterSetRNGStream(
    cl,
    iseed = seed + 1
  )
  
  doParallel::registerDoParallel(
    cl
  )
  
  
  ##########################################################
  ## Parallel simulation
  ##########################################################
  
  results <- foreach(
    
    k =
      seq_len(
        nrow(grid)
      ),
    
    .combine =
      rbind,
    
    .packages =
      character(0),
    
    .export =
      c(
        "legendre_poly",
        "phi_basis",
        
        "generate_fixed_theta",
        "compute_f_matrix",
        "get_Bprime",
        
        "solve_constrained_from_eigen",
        "prepare_design",
        
        "generate_noisy_dataset",
        
        "prepare_A1",
        "privatize_A1",
        
        "prepare_A3",
        "privatize_A3",
        
        "prepare_PTF",
        "privatize_PTF",
        
        "prepare_PTF_WLS",
        "fit_A2",
        "fit_A4",
        
        "compute_ISE",
        
        "compute_IVar_FTP",
        "compute_IVar_PTF",
        
        "run_one_setting"
      )
  ) %dopar% {
    
    ############################################################
    ## Current setting
    ############################################################
    
    n <- grid$n[k]
    m_current <- grid$m[k]
    T_current <- grid$T[k]
    B_current <- grid$B[k]
    sigma_current <- grid$sigma_obs[k]
    
    Bprime_choice <- grid$Bprime_choice[k]
    mu <- grid$mu[k]
    
    ############################################################
    ## Retrieve fixed truth and design
    ############################################################
    
    truth_key <- paste(
      n, m_current, B_current, sep = "_"
    )
    
    design_key <- paste(
      T_current, m_current, sep = "_"
    )
    
    truth <- fixed_data[[truth_key]]
    
    current_design <- design_data[[design_key]]
    
    Phi <- current_design$Phi
    design <- current_design$design
    phi_norm <- current_design$phi_norm
    PTF_WLS <- current_design$PTF_WLS
    
    ############################################################
    ## Evaluate fixed latent functions on this grid
    ############################################################
    
    f_mat <- compute_f_matrix(
      truth$theta,
      Phi
    )
    
    ############################################################
    ## Calculate B'
    ############################################################
    
    Bprime <- get_Bprime(
      n = n,
      B = B_current,
      choice = Bprime_choice
    )
    
    ############################################################
    ## Run Monte Carlo simulation
    ############################################################
    
    run_one_setting(
      theta = truth$theta,
      f_mat = f_mat,
      theta_bar_true = truth$theta_bar_true,
      
      Phi = Phi,
      design = design,
      PTF_WLS = PTF_WLS,
      
      n = n,
      T = T_current,
      m = m_current,
      
      mu = mu,
      B = B_current,
      Bprime = Bprime,
      Bprime_choice = Bprime_choice,
      
      sigma_obs = sigma_current,
      phi_norm = phi_norm,
      
      reps = reps
    )
  }
  
  
  ##########################################################
  ## Clean row names
  ##########################################################
  
  rownames(results) <- NULL
  
  cat(
    "Simulation complete.\n"
  )
  
  results
}


############################################################
##### Run simulation
############################################################

tbl_all <- run_simulation(
  n_vec = 2^(2:11),
  mu_vec = c(1),
  m_vec = c(2, 3, 4, 5, 6, 8, 10),
  Bprime_choices = c("sqrtlog"),
  T = c(20),
  B = c(1),
  mu0 = NULL,
  sigma_obs = c(0),
  reps = 5000,
  seed = 12345
)

##### Simulation 1 for Noiseless

#tbl_all <- run_simulation(
#  n_vec = c(5, 6, 7, 8, 9, 10, 15, 20, 30), mu = 1
#  n_vec = c(50, 60, 70, 80, 90, 100, 150, 200, 300), mu = 0.1
#  m_vec = c(2, 3, 4, 5, 6, 8, 10),
#  Bprime_choices = c("log"),
#  T = c(20),
#  B = c(1),
#  mu0 = NULL,
#  sigma_obs = c(0),
#  reps = 5000,
#  seed = 12345
#)

##### Simulation 2 for B' choices
# tbl_all <- run_simulation(
#   n_vec = c(4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048),
#   mu_vec = c(1),
#   m_vec = c(4),
#   Bprime_choices = c("B", "sqrtlog", "log", "sqrtn", "n_2_3"),
#   T = c(20),
#   B = c(1),
#   mu0 = NULL,
#   sigma_obs = c(1),
#   reps = 5000,
#   seed = 12345
# )

##### Simulation 3 for Noisy Level
#tbl_all <- run_simulation(
#  n_vec = c(4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048),
#  mu_vec = c(1),
#  m_vec = c(2, 3, 4, 5, 6, 8, 10),
#  Bprime_choices = c("log"),
#  T = c(20),
#  B = c(1),
#  mu0 = NULL,
#  sigma_obs = c(0.01, 1, 100),
#  reps = 5000,
#  seed = 12345
#)

############################################################
##### Sort results
############################################################

tbl_all <-
  tbl_all[
    order(
      tbl_all$n,
      tbl_all$m,
      tbl_all$T,
      tbl_all$B,
      tbl_all$sigma_obs,
      tbl_all$mu,
      match(
        tbl_all$Bprime_choice,
        c(
          "B",
          "sqrtlog",
          "log",
          "sqrtn",
          "n_2_3"
        )
      )
    ),
  ]

rownames(tbl_all) <- NULL

############################################################
##### Save
############################################################

saveRDS(
  tbl_all,
  file =
    "Poly_sigma0_0_mu_1_n_4-2048.rds"
)
