# ================================================================
#  Section 1: Simulation Setup
#
#  Scenarios:
#    S_base : Smooth Matérn coefficient surfaces
#    T_d05  : Soft boundary transition, delta = 0.5
#    T_d10  : Soft boundary transition, delta = 1.0
#    T_d15  : Soft boundary transition, delta = 1.5
#    H_d05  : Hard step discontinuity, delta = 0.5
#    H_d10  : Hard step discontinuity, delta = 1.0
#    H_d15  : Hard step discontinuity, delta = 1.5
#    H_d25  : Hard step discontinuity, delta = 2.5
#
#  Output:
#    sim_fixed_master.RData
# ================================================================

Sys.setlocale("LC_ALL", "C")

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(MASS)
})

# ================================================================
# PART 0: Global configuration
# ================================================================

N_SIDE    <- 25L
N_LOC     <- N_SIDE^2L
M_REP     <- 3L
N_OBS     <- N_LOC * M_REP
NB_SIZE   <- 8L

BASE_SEED <- 20260301L

# Monte Carlo design
N_REP_MC <- 200L

# Scenario intensity levels
DELTA_SOFT <- c(d05 = 0.5, d10 = 1.0, d15 = 1.5)
DELTA_HARD <- c(d05 = 0.5, d10 = 1.0, d15 = 1.5, d25 = 2.5)

# Soft transition width
TAU_SOFT <- 1.0

# ================================================================
# PART 1: Smooth Matérn background parameters
# ================================================================

# beta0(s): intercept
RANGE0_BG <- 8.0
SIGMA0_BG <- 1.0
MEAN0_BG  <- 2.0

# beta1(s): slope for x1
RANGE1_BG <- 12.0
SIGMA1_BG <- 0.6
MEAN1_BG  <- 0.5

# beta2(s): slope for x2, smooth control coefficient
RANGE2_BG <- 6.0
SIGMA2_BG <- 0.5
MEAN2_BG  <- -0.3

# ================================================================
# PART 2: MGWNBR prior and kNN configuration
# ================================================================

# Maximum allowed neighbourhood size
K_MAX <- as.integer(N_LOC * 0.40)

# Precomputed neighbour order
K_PRE <- min(120L, N_LOC - 1L)

# Coefficient-specific priors for k
B0_KSHAPE <- 4.0
K0_CTR    <- 50L
K_MIN_B0  <- 15L

B1_KSHAPE <- 4.0
K1_CTR    <- 60L
K_MIN_B1  <- 20L

B2_KSHAPE <- 5.0
K2_CTR    <- 40L
K_MIN_B2  <- 12L

K0_RATE <- K0_CTR * (B0_KSHAPE + 1)
K1_RATE <- K1_CTR * (B1_KSHAPE + 1)
K2_RATE <- K2_CTR * (B2_KSHAPE + 1)

# PC-type prior on latent marginal SD through tau
V0 <- 1.0
V1 <- 0.5
V2 <- 0.8

ALPHA_PC_B0 <- 0.5
ALPHA_PC_B1 <- 0.5
ALPHA_PC_B2 <- 0.5

LAMBDA_B0 <- -log(ALPHA_PC_B0) / sqrt(V0)
LAMBDA_B1 <- -log(ALPHA_PC_B1) / sqrt(V1)
LAMBDA_B2 <- -log(ALPHA_PC_B2) / sqrt(V2)

TH1_B0 <- 2 * log(LAMBDA_B0)
TH1_B1 <- 2 * log(LAMBDA_B1)
TH1_B2 <- 2 * log(LAMBDA_B2)

# ================================================================
# PART 3: Spatial coordinates
# ================================================================

uv <- expand.grid(
  east  = seq(0, 24, length.out = N_SIDE),
  north = seq(0, 24, length.out = N_SIDE)
)

east   <- uv$east
north  <- uv$north
coords <- cbind(east, north)

cx <- 12.0
cy <- 12.0

idx_rep <- rep(seq_len(N_LOC), each = M_REP)
x0_rep  <- rep(1.0, N_OBS)

Dist0 <- as.matrix(dist(coords))
DMAX  <- max(Dist0)

D <- Dist0 / DMAX

cat(sprintf("[Coordinates] N_LOC=%d | DMAX=%.4f\n\n", N_LOC, DMAX))

# ================================================================
# PART 4: generate Matérn-3/2 Gaussian fields
# ================================================================

gen_matern_gp <- function(coords, range_val, sigma_val, mean_val = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n   <- nrow(coords)
  D_m <- as.matrix(dist(coords))
  
  s3  <- sqrt(3) * D_m / range_val
  Sig <- sigma_val^2 * (1 + s3) * exp(-s3)
  
  # Small jitter for numerical stability
  diag(Sig) <- diag(Sig) + 1e-6
  
  mean_val + as.vector(t(chol(Sig)) %*% rnorm(n))
}

# ================================================================
# PART 5: Matérn background fields
# ================================================================

cat("[GP] Generating fixed Matérn background fields...\n")

beta0_smooth <- gen_matern_gp(
  coords     = coords,
  range_val  = RANGE0_BG,
  sigma_val  = SIGMA0_BG,
  mean_val   = MEAN0_BG,
  seed       = 201
)

beta1_smooth <- gen_matern_gp(
  coords     = coords,
  range_val  = RANGE1_BG,
  sigma_val  = SIGMA1_BG,
  mean_val   = MEAN1_BG,
  seed       = 202
)

beta2_smooth <- gen_matern_gp(
  coords     = coords,
  range_val  = RANGE2_BG,
  sigma_val  = SIGMA2_BG,
  mean_val   = MEAN2_BG,
  seed       = 203
)

beta2_true <- beta2_smooth

# ================================================================
# PART 6: Boundary geometry
# ================================================================

# ------------------------------------------------
# Boundary for beta0:
# vertical-horizontal cross centred at (12,12)
# ------------------------------------------------
z0_hard <- ifelse(
  east <= cx & north <= cy, -1.0,
  ifelse(
    east <= cx & north > cy, +1.0,
    ifelse(east > cx & north <= cy, +1.0, -1.0)
  )
)

# Distance to nearest vertical or horizontal boundary
dist_to_vh <- pmin(abs(east - cx), abs(north - cy))
d0_signed <- z0_hard * dist_to_vh

# ------------------------------------------------
# Boundary for beta1:
# anti-diagonal e + n = 24
# ------------------------------------------------

diag_flag <- (east + north) < 24.0

z1_hard <- ifelse(diag_flag, +1.0, -1.0)

dist_to_diag <- abs(east + north - 24.0) / sqrt(2)

d1_signed <- z1_hard * dist_to_diag

# ------------------------------------------------
# Standardised hard step functions
# ------------------------------------------------

jump0_hard_raw <- z0_hard
jump1_hard_raw <- z1_hard

jump0_hard <- jump0_hard_raw / sd(jump0_hard_raw)
jump1_hard <- jump1_hard_raw / sd(jump1_hard_raw)

# ------------------------------------------------
# Soft transition functions
# ------------------------------------------------

jump0_soft_raw <- tanh(d0_signed / TAU_SOFT)
jump1_soft_raw <- tanh(d1_signed / TAU_SOFT)

jump0_soft <- jump0_soft_raw / sd(jump0_soft_raw)
jump1_soft <- jump1_soft_raw / sd(jump1_soft_raw)

# ================================================================
# PART 7: Scenario definitions and true coefficient surfaces
# ================================================================

scenario_table <- data.frame(
  scenario_id = c(
    "S_base",
    paste0("T_", names(DELTA_SOFT)),
    paste0("H_", names(DELTA_HARD))
  ),
  scenario_family = c(
    "S",
    rep("T", length(DELTA_SOFT)),
    rep("H", length(DELTA_HARD))
  ),
  scenario_label = c(
    "Smooth Matern",
    rep("Soft transition", length(DELTA_SOFT)),
    rep("Hard discontinuity", length(DELTA_HARD))
  ),
  delta = c(
    NA_real_,
    as.numeric(DELTA_SOFT),
    as.numeric(DELTA_HARD)
  ),
  delta_name = c(
    "base",
    names(DELTA_SOFT),
    names(DELTA_HARD)
  ),
  stringsAsFactors = FALSE
)

beta_true_list <- list()

# ------------------------------------------------
# Scenario S: Smooth baseline
# ------------------------------------------------

beta_true_list[["S_base"]] <- list(
  beta0 = beta0_smooth,
  beta1 = beta1_smooth,
  beta2 = beta2_true,
  scenario_family = "S",
  scenario_label  = "Smooth Matern",
  delta           = NA_real_,
  delta_name      = "base",
  jump0           = rep(0, N_LOC),
  jump1           = rep(0, N_LOC)
)

# ------------------------------------------------
# Scenario T: Soft transitions
# ------------------------------------------------

for (nm in names(DELTA_SOFT)) {
  delta_val <- DELTA_SOFT[[nm]]
  sid <- paste0("T_", nm)
  
  beta_true_list[[sid]] <- list(
    beta0 = beta0_smooth + delta_val * jump0_soft,
    beta1 = beta1_smooth + delta_val * jump1_soft,
    beta2 = beta2_true,
    scenario_family = "T",
    scenario_label  = "Soft transition",
    delta           = delta_val,
    delta_name      = nm,
    jump0           = jump0_soft,
    jump1           = jump1_soft
  )
}

# ------------------------------------------------
# Scenario H: Hard discontinuities
# ------------------------------------------------

for (nm in names(DELTA_HARD)) {
  delta_val <- DELTA_HARD[[nm]]
  sid <- paste0("H_", nm)
  
  beta_true_list[[sid]] <- list(
    beta0 = beta0_smooth + delta_val * jump0_hard,
    beta1 = beta1_smooth + delta_val * jump1_hard,
    beta2 = beta2_true,
    scenario_family = "H",
    scenario_label  = "Hard discontinuity",
    delta           = delta_val,
    delta_name      = nm,
    jump0           = jump0_hard,
    jump1           = jump1_hard
  )
}

# ================================================================
# PART 8: Boundary masks and distance bins
# ================================================================

# Near-boundary masks used for simple near/far diagnostics
near_vhbnd <- abs(east - cx) < 2 | abs(north - cy) < 2
near_diag  <- abs(east + north - 24.0) < 3

# Distance bins for boundary-distance RMSE decomposition
BND_BREAKS <- c(0, 1, 2, 4, 8, Inf)
BND_LABELS <- c("[0,1)", "[1,2)", "[2,4)", "[4,8)", "[8,Inf)")

bnd_bin_vh <- cut(
  dist_to_vh,
  breaks = BND_BREAKS,
  labels = BND_LABELS,
  right = FALSE,
  include.lowest = TRUE
)

bnd_bin_diag <- cut(
  dist_to_diag,
  breaks = BND_BREAKS,
  labels = BND_LABELS,
  right = FALSE,
  include.lowest = TRUE
)

# ================================================================
# PART 9: kNN precomputation for MGWNBR
# ================================================================

t_pre <- proc.time()

knn_idx_mat <- matrix(0L,  N_LOC, K_PRE)
D_knn_mat   <- matrix(0.0, N_LOC, K_PRE)

for (i in seq_len(N_LOC)) {
  ord_i <- order(D[i, ])[2:(K_PRE + 1)]
  knn_idx_mat[i, ] <- ord_i
  D_knn_mat[i, ]   <- D[i, ord_i]
}

sorted_D_mat <- apply(D, 1, function(row) sort(row[row > 1e-10]))

med_bw_at_k_vec <- sapply(seq_len(K_PRE), function(ki) {
  median(sorted_D_mat[ki, ]) * DMAX
})

# ================================================================
# PART 10: SPDE mesh and objects for SVC-NBR
# ================================================================

cat("[SPDE] Building mesh...\n")

mesh_sim <- inla.mesh.2d(
  loc      = coords,
  max.edge = c(2.0, 6.0),
  offset   = c(1.5, 4.0),
  cutoff   = 1.0
)

cat(sprintf("[SPDE] mesh nodes = %d\n", mesh_sim$n))

# PC priors
spde_sim_b0 <- inla.spde2.pcmatern(
  mesh_sim,
  alpha       = 2,
  prior.range = c(4.0, 0.05),
  prior.sigma = c(2.0, 0.05)
)

spde_sim_b1 <- inla.spde2.pcmatern(
  mesh_sim,
  alpha       = 2,
  prior.range = c(4.0, 0.05),
  prior.sigma = c(1.5, 0.05)
)

spde_sim_b2 <- inla.spde2.pcmatern(
  mesh_sim,
  alpha       = 2,
  prior.range = c(3.0, 0.05),
  prior.sigma = c(1.5, 0.05)
)

A_sim     <- inla.spde.make.A(mesh_sim, loc = coords)
A_sim_rep <- A_sim[idx_rep, ]

idx0_svc <- inla.spde.make.index("s0", n.spde = spde_sim_b0$n.spde)
idx1_svc <- inla.spde.make.index("s1", n.spde = spde_sim_b1$n.spde)
idx2_svc <- inla.spde.make.index("s2", n.spde = spde_sim_b2$n.spde)

cat("[SPDE] Objects created.\n\n")

# ================================================================
# PART 11: Save setup
# ================================================================

save(
  # Global simulation settings
  N_SIDE, N_LOC, M_REP, N_OBS, NB_SIZE,
  BASE_SEED, N_REP_MC,
  DELTA_SOFT, DELTA_HARD, TAU_SOFT,
  
  # Background parameters
  RANGE0_BG, SIGMA0_BG, MEAN0_BG,
  RANGE1_BG, SIGMA1_BG, MEAN1_BG,
  RANGE2_BG, SIGMA2_BG, MEAN2_BG,
  
  # MGWNBR priors and k configuration
  K_MAX, K_PRE,
  K0_CTR, K_MIN_B0, B0_KSHAPE, K0_RATE, LAMBDA_B0, TH1_B0,
  K1_CTR, K_MIN_B1, B1_KSHAPE, K1_RATE, LAMBDA_B1, TH1_B1,
  K2_CTR, K_MIN_B2, B2_KSHAPE, K2_RATE, LAMBDA_B2, TH1_B2,
  V0, V1, V2,
  ALPHA_PC_B0, ALPHA_PC_B1, ALPHA_PC_B2,
  
  # Coordinates and distances
  east, north, coords, cx, cy,
  Dist0, DMAX, D,
  idx_rep, x0_rep,
  
  # Smooth background fields
  beta0_smooth, beta1_smooth, beta2_smooth, beta2_true,
  
  # Boundary objects
  diag_flag,
  z0_hard, z1_hard,
  dist_to_vh, dist_to_diag,
  d0_signed, d1_signed,
  jump0_hard_raw, jump1_hard_raw,
  jump0_hard, jump1_hard,
  jump0_soft_raw, jump1_soft_raw,
  jump0_soft, jump1_soft,
  near_vhbnd, near_diag,
  BND_BREAKS, BND_LABELS,
  bnd_bin_vh, bnd_bin_diag,
  
  # Scenarios and true coefficients
  scenario_table, beta_true_list,
  
  # kNN precomputation
  knn_idx_mat, D_knn_mat, sorted_D_mat, med_bw_at_k_vec,
  
  # SPDE objects
  mesh_sim,
  spde_sim_b0, spde_sim_b1, spde_sim_b2,
  A_sim, A_sim_rep,
  idx0_svc, idx1_svc, idx2_svc,
  
  file = "sim_fixed_master.RData"
)


# ================================================================
#  Section 2: run the simulations via HPC
#
#  Arguments:
#    kk          : Monte Carlo replicate index, e.g. 1,...,200
#    scenario_id : one of
#                  S_base,
#                  T_d05, T_d10, T_d15,
#                  H_d05, H_d10, H_d15, H_d25
#
#  Required input:
#    sim_fixed_master.RData
#
#  Output:
#    results_master/<scenario_id>/<kk>.RData
# ================================================================

Sys.setlocale("LC_ALL", "C")

if (dir.exists("~/rlibs")) {
  .libPaths("~/rlibs")
}

args <- commandArgs(trailingOnly = TRUE)

kk          <- as.integer(args[1])
scenario_id <- as.character(args[2])
force_rerun <- as.integer(args[3])


suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(MASS)
})

# INLA options
inla.setOption(num.threads = "4:1")
try(inla.setOption(smtp = "taucs"), silent = TRUE)

# ================================================================
# PART 1: Load master setup
# ================================================================

setup_file <- "sim_fixed_master.RData"
load(setup_file)

valid_scenarios <- scenario_table$scenario_id

if (!scenario_id %in% valid_scenarios) {
  stop(sprintf(
    "Invalid scenario_id: %s. Valid values are: %s",
    scenario_id,
    paste(valid_scenarios, collapse = ", ")
  ))
}

scenario_info <- scenario_table[scenario_table$scenario_id == scenario_id, ]
beta_true_obj <- beta_true_list[[scenario_id]]

beta0_true <- beta_true_obj$beta0
beta1_true <- beta_true_obj$beta1
beta2_true <- beta_true_obj$beta2

# ================================================================
# PART 2: Output path and skip rule
# ================================================================

out_root <- "results_master"
out_dir  <- file.path(out_root, scenario_id)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

out_file <- file.path(out_dir, sprintf("%d.RData", kk))

# ================================================================
# PART 3: Helper functions
# ================================================================

mk_constr <- function(sp) {
  list(A = matrix(1, 1, sp$n.spde), e = 0)
}

mk_constr_n <- function(nn) {
  list(A = matrix(1, 1, nn), e = 0)
}

safe_extract <- function(x, nm, default = NA_real_) {
  if (is.null(x)) return(default)
  if (!nm %in% names(x)) return(default)
  x[[nm]]
}

calc_metrics <- function(est, true, lo, hi) {
  err <- est - true
  
  bias   <- mean(err, na.rm = TRUE)
  sd_err <- sd(err, na.rm = TRUE)
  rmse   <- sqrt(mean(err^2, na.rm = TRUE))
  mae    <- mean(abs(err), na.rm = TRUE)
  
  r2 <- if (var(true, na.rm = TRUE) < 1e-10) {
    NA_real_
  } else {
    1 - sum((est - true)^2, na.rm = TRUE) /
      sum((true - mean(true, na.rm = TRUE))^2, na.rm = TRUE)
  }
  
  cov95 <- mean(true >= lo & true <= hi, na.rm = TRUE)
  ci_w  <- mean(hi - lo, na.rm = TRUE)
  
  calib <- if (!is.finite(rmse) || rmse <= 1e-12) {
    NA_real_
  } else {
    mean((hi - lo) / (2 * 1.96), na.rm = TRUE) / rmse
  }
  
  list(
    bias        = bias,
    sd_err      = sd_err,
    rmse        = rmse,
    mae         = mae,
    r2          = r2,
    cov95       = cov95,
    ci_width    = ci_w,
    calib_ratio = calib
  )
}

bnd_rmse <- function(err_vec, mask) {
  near_rmse <- sqrt(mean(err_vec[mask]^2, na.rm = TRUE))
  far_rmse  <- sqrt(mean(err_vec[!mask]^2, na.rm = TRUE))
  
  list(
    near  = near_rmse,
    far   = far_rmse,
    ratio = near_rmse / max(far_rmse, 1e-8)
  )
}

bnd_profile <- function(err_vec,
                        dist_vec,
                        breaks = BND_BREAKS,
                        labels = BND_LABELS) {
  cuts <- cut(
    dist_vec,
    breaks = breaks,
    labels = labels,
    right = FALSE,
    include.lowest = TRUE
  )
  
  out <- tapply(err_vec^2, cuts, function(x) {
    sqrt(mean(x, na.rm = TRUE))
  })
  
  # Ensure all bins are present in a fixed order.
  out_full <- setNames(rep(NA_real_, length(labels)), labels)
  out_full[names(out)] <- as.numeric(out)
  
  out_full
}

extract_spde_hp <- function(hp, fname) {
  if (is.null(hp) || nrow(hp) == 0L) {
    return(list(
      range_mean = NA_real_,
      range_q025 = NA_real_,
      range_q975 = NA_real_,
      sigma_mean = NA_real_,
      sigma_q025 = NA_real_,
      sigma_q975 = NA_real_
    ))
  }
  
  rn <- rownames(hp)
  
  r_row <- grep(sprintf("Range.*%s", fname), rn, ignore.case = TRUE)
  s_row <- grep(sprintf("Stdev.*%s", fname), rn, ignore.case = TRUE)
  
  out <- list(
    range_mean = NA_real_,
    range_q025 = NA_real_,
    range_q975 = NA_real_,
    sigma_mean = NA_real_,
    sigma_q025 = NA_real_,
    sigma_q975 = NA_real_
  )
  
  if (length(r_row)) {
    rr <- r_row[1]
    out$range_mean <- hp[rr, "mean"]
    out$range_q025 <- hp[rr, "0.025quant"]
    out$range_q975 <- hp[rr, "0.975quant"]
  }
  
  if (length(s_row)) {
    ss <- s_row[1]
    out$sigma_mean <- hp[ss, "mean"]
    out$sigma_q025 <- hp[ss, "0.025quant"]
    out$sigma_q975 <- hp[ss, "0.975quant"]
  }
  
  out
}

get_theta2 <- function(hp, tag) {
  if (is.null(hp) || nrow(hp) == 0L) {
    return(list(
      k_mean  = NA_real_,
      k_q025  = NA_real_,
      k_q975  = NA_real_,
      bw_mean = NA_real_
    ))
  }
  
  row <- grep(sprintf("Theta2.*%s", tag), rownames(hp), ignore.case = TRUE)
  
  if (!length(row)) {
    return(list(
      k_mean  = NA_real_,
      k_q025  = NA_real_,
      k_q975  = NA_real_,
      bw_mean = NA_real_
    ))
  }
  
  rr <- row[1]
  
  k_mean <- exp(hp[rr, "mean"])
  k_q025 <- exp(hp[rr, "0.025quant"])
  k_q975 <- exp(hp[rr, "0.975quant"])
  
  ki <- max(1L, min(round(k_mean), length(med_bw_at_k_vec)))
  
  list(
    k_mean  = k_mean,
    k_q025  = k_q025,
    k_q975  = k_q975,
    bw_mean = med_bw_at_k_vec[ki]
  )
}

get_lcpo <- function(fit) {
  if (is.null(fit$cpo) || is.null(fit$cpo$cpo)) {
    return(NA_real_)
  }
  
  cpo_v <- fit$cpo$cpo
  
  valid <- cpo_v > 0 &
    is.finite(cpo_v) &
    !is.na(cpo_v)
  
  if (!sum(valid)) {
    NA_real_
  } else {
    -mean(log(cpo_v[valid]))
  }
}

get_nb_size <- function(hp) {
  if (is.null(hp) || nrow(hp) == 0L) {
    return(NA_real_)
  }
  
  nb_row <- grep("size|nbinom|overdispersion", rownames(hp), ignore.case = TRUE)
  
  if (!length(nb_row)) {
    return(NA_real_)
  }
  
  hp[nb_row[1], "mean"]
}

get_ic <- function(fit) {
  list(
    dic  = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
    waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_,
    lcpo = get_lcpo(fit),
    mlik = if (!is.null(fit$mlik)) as.numeric(fit$mlik[1, 1]) else NA_real_
  )
}

safe_fixed_mean <- function(fe, nm) {
  if (is.null(fe)) return(0)
  if (!nm %in% rownames(fe)) return(0)
  fe[nm, "mean"]
}

# ================================================================
# PART 4: MGWNBR precision model
# ================================================================

inla.rgeneric.mgwnbr.master <- function(
    cmd   = c("graph", "Q", "mu", "initial",
              "log.norm.const", "log.prior", "quit"),
    theta = NULL) {
  
  envir <- parent.env(environment())
  
  itp <- function() {
    tau <- exp(max(min(theta[1L], 12.0), -4.0))
    
    k_raw <- exp(
      max(
        min(theta[2L], log(K_MAX)),
        log(K_MIN)
      )
    )
    
    list(tau = tau, k = k_raw)
  }
  
  graph <- function() {
    K_G <- min(K_PRE_G, 60L)
    
    i_idx <- rep(seq_len(n_g), each = K_G)
    j_idx <- as.vector(t(KNN_IDX_G[, seq_len(K_G)]))
    
    all_i <- c(i_idx, j_idx, seq_len(n_g))
    all_j <- c(j_idx, i_idx, seq_len(n_g))
    
    G <- Matrix::sparseMatrix(
      i = all_i,
      j = all_j,
      x = 1.0,
      dims = c(n_g, n_g),
      use.last.ij = TRUE
    )
    
    as(as(G, "generalMatrix"), "dgCMatrix")
  }
  
  Q <- function() {
    p <- itp()
    
    k_lo <- max(1L, floor(p$k))
    k_hi <- min(K_PRE_G, ceiling(p$k))
    al   <- p$k - k_lo
    
    bw_lo <- KNN_D_G[k_lo, ]
    bw_hi <- KNN_D_G[k_hi, ]
    
    b_vec <- (1.0 - al) * bw_lo + al * bw_hi
    b_vec <- pmax(b_vec, 1e-4)
    
    K_USE <- min(K_PRE_G, max(ceiling(p$k) + 8L, 25L))
    
    i_idx  <- rep(seq_len(n_g), each = K_USE)
    j_idx  <- as.vector(t(KNN_IDX_G[, seq_len(K_USE)]))
    d_vals <- as.vector(t(KNN_D_G[seq_len(K_USE), ]))
    b_src  <- rep(b_vec, each = K_USE)
    
    w_vals <- exp(-d_vals^2 / (2.0 * b_src^2))
    w_vals[!is.finite(w_vals)] <- 0.0
    w_vals[w_vals < 1e-7]      <- 0.0
    
    W_sp <- Matrix::sparseMatrix(
      i = i_idx,
      j = j_idx,
      x = w_vals,
      dims = c(n_g, n_g),
      use.last.ij = FALSE
    )
    
    W_sym <- (W_sp + Matrix::t(W_sp)) / 2.0
    
    rs <- Matrix::rowSums(W_sym)
    
    # Small stabilising nugget.
    nug <- p$tau * max(rs, 1.0) * 1e-2 + 1e-6
    
    Q_out <- p$tau * (Matrix::Diagonal(x = rs) - W_sym) +
      nug * Matrix::Diagonal(n_g)
    
    as(as(Q_out, "CsparseMatrix"), "dgCMatrix")
  }
  
  mu <- function() {
    rep(0.0, n_g)
  }
  
  log.norm.const <- function() {
    numeric(0)
  }
  
  log.prior <- function() {
    tau <- exp(theta[1L])
    
    if (!is.finite(tau) || tau <= 0) {
      return(-1e10)
    }
    
    sigma <- 1.0 / sqrt(tau)
    
    # PC-type prior on sigma via tau.
    lp_tau <- log(LAMBDA_G / 2) -
      theta[1L] / 2 -
      LAMBDA_G * sigma
    
    # Prior on k through theta2 = log(k).
    lp_k <- -(K_SHAPE_G + 1) * theta[2L] -
      K_RATE_G * exp(-theta[2L])
    
    lp <- lp_tau + lp_k
    
    if (!is.finite(lp)) {
      return(-1e10)
    }
    
    lp
  }
  
  initial <- function() {
    c(TH1_INIT_G, log(K_CTR_G))
  }
  
  quit <- function() {
    invisible()
  }
  
  if (!length(theta) || is.null(theta)) {
    theta <- initial()
  }
  
  do.call(match.arg(cmd), args = list())
}

mk_model_mg <- function(th1, k_ctr, k_min, k_shape, k_rate, lam) {
  inla.rgeneric.define(
    inla.rgeneric.mgwnbr.master,
    n_g        = N_LOC,
    KNN_IDX_G  = knn_idx_mat,
    KNN_D_G    = t(D_knn_mat),
    K_PRE_G    = K_PRE,
    TH1_INIT_G = th1,
    K_CTR_G    = k_ctr,
    K_MIN      = k_min,
    K_MAX      = K_MAX,
    K_SHAPE_G  = k_shape,
    K_RATE_G   = k_rate,
    LAMBDA_G   = lam
  )
}

# ================================================================
# PART 5: MGWNBR fitting
# ================================================================

fit_mgwnbr_robust <- function(f, data_df, nb_prior, verbose = FALSE) {
  
  strategies <- list(
    list(
      tag = "Gaussian+EB",
      ci  = list(
        strategy     = "gaussian",
        int.strategy = "eb",
        diagonal     = 1.0,
        h            = 0.01,
        tolerance    = 1e-3
      ),
      full = FALSE
    ),
    list(
      tag = "SimpLaplace+EB",
      ci  = list(
        strategy     = "simplified.laplace",
        int.strategy = "eb",
        diagonal     = 1e-2,
        h            = 8e-3,
        tolerance    = 1e-4
      ),
      full = TRUE
    ),
    list(
      tag = "SimpLaplace+EB_safe",
      ci  = list(
        strategy     = "simplified.laplace",
        int.strategy = "eb",
        diagonal     = 0.1,
        h            = 0.02,
        tolerance    = 1e-3
      ),
      full = TRUE
    )
  )
  
  theta_mode <- NULL
  m_final    <- NULL
  fit_log    <- list()
  
  for (s in strategies) {
    cat(sprintf("  -> MGWNBR strategy: %s\n", s$tag))
    
    ctrl_mode <- if (!is.null(theta_mode)) {
      list(theta = theta_mode, restart = TRUE)
    } else {
      list(result = NULL)
    }
    
    m_try <- tryCatch(
      inla(
        f,
        family = "nbinomial",
        data = data_df,
        safe = TRUE,
        control.family = nb_prior,
        control.predictor = list(compute = s$full),
        control.compute = if (s$full) {
          list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE)
        } else {
          list(dic = FALSE, waic = FALSE, cpo = FALSE, mlik = FALSE)
        },
        control.mode = ctrl_mode,
        control.inla = s$ci,
        verbose = verbose
      ),
      error = function(e) {
        msg <- conditionMessage(e)
        cat(sprintf("     FAIL: %s\n", substr(msg, 1, 120)))
        fit_log[[s$tag]] <<- list(status = "FAIL", message = msg)
        NULL
      }
    )
    
    if (!is.null(m_try)) {
      theta_mode <- m_try$mode$theta
      
      cat(sprintf(
        "     OK | theta=%s\n",
        paste(round(theta_mode, 3), collapse = ", ")
      ))
      
      fit_log[[s$tag]] <- list(
        status = "OK",
        theta  = theta_mode
      )
      
      if (s$full) {
        m_final <- m_try
        break
      }
    }
  }
  
  list(fit = m_final, log = fit_log)
}

# ================================================================
# PART 6: Data generation
# ================================================================

seed_this <- BASE_SEED + kk

# Scenario-specific offset in seed to avoid identical covariates/responses
# across scenarios if desired. This preserves reproducibility.
scenario_offset <- match(scenario_id, valid_scenarios) * 100000L
seed_data <- seed_this + scenario_offset

set.seed(seed_data)

cat("[Data generation]\n")
cat(sprintf("  seed_data = %d\n", seed_data))

x1_rep <- rnorm(N_OBS, 0, 1)
x2_rep <- rnorm(N_OBS, 0, 1)

b0_rep <- rep(beta0_true, each = M_REP)
b1_rep <- rep(beta1_true, each = M_REP)
b2_rep <- rep(beta2_true, each = M_REP)

eta_true <- b0_rep +
  b1_rep * x1_rep +
  b2_rep * x2_rep

# Avoid numerical overflow in extreme stress tests.
eta_true_clip <- pmin(pmax(eta_true, -20), 20)

mu_true <- exp(eta_true_clip)

y_rep <- rnbinom(
  n    = N_OBS,
  size = NB_SIZE,
  mu   = mu_true
)

count_diag <- summarise_counts(y_rep)

# ================================================================
# PART 7: Result object initialisation
# ================================================================

RSLT <- list(
  kk          = kk,
  seed_data   = seed_data,
  scenario_id = scenario_id,
  scenario    = as.list(scenario_info),
  count_diag  = count_diag,
  
  mgwnbr = NULL,
  svc    = NULL,
  
  mgwnbr_error = NULL,
  svc_error    = NULL,
  
  fit_status = list(
    mgwnbr_ok = FALSE,
    svc_ok    = FALSE,
    both_ok   = FALSE
  ),
  
  comparison = NULL
)

# ================================================================
# PART 8: MODEL 1 -- MGWNBR
# ================================================================

t0_mg <- proc.time()

m_mg <- NULL
mg_fit_log <- NULL

tryCatch({
  
  model.b0 <- mk_model_mg(
    TH1_B0,
    K0_CTR,
    K_MIN_B0,
    B0_KSHAPE,
    K0_RATE,
    LAMBDA_B0
  )
  
  model.b1 <- mk_model_mg(
    TH1_B1,
    K1_CTR,
    K_MIN_B1,
    B1_KSHAPE,
    K1_RATE,
    LAMBDA_B1
  )
  
  model.b2 <- mk_model_mg(
    TH1_B2,
    K2_CTR,
    K_MIN_B2,
    B2_KSHAPE,
    K2_RATE,
    LAMBDA_B2
  )
  
  data_mg <- data.frame(
    y    = y_rep,
    x0   = x0_rep,
    x1   = x1_rep,
    x2   = x2_rep,
    
    idx0 = idx_rep,
    idx1 = idx_rep,
    idx2 = idx_rep,
    
    mu0  = x0_rep,
    mu1  = x1_rep,
    mu2  = x2_rep
  )
  
  f_mg <- y ~ -1 + mu0 + mu1 + mu2 +
    f(
      idx0,
      x0,
      model = model.b0,
      n = N_LOC,
      extraconstr = mk_constr_n(N_LOC)
    ) +
    f(
      idx1,
      x1,
      model = model.b1,
      n = N_LOC,
      extraconstr = mk_constr_n(N_LOC)
    ) +
    f(
      idx2,
      x2,
      model = model.b2,
      n = N_LOC,
      extraconstr = mk_constr_n(N_LOC)
    )
  
  nb_prior_mg <- list(
    hyper = list(
      theta = list(
        prior   = "normal",
        param   = c(log(NB_SIZE), 0.15),
        initial = log(NB_SIZE)
      )
    )
  )
  
  mg_fit_obj <- fit_mgwnbr_robust(
    f        = f_mg,
    data_df  = data_mg,
    nb_prior = nb_prior_mg,
    verbose  = FALSE
  )
  
  m_mg       <- mg_fit_obj$fit
  mg_fit_log <- mg_fit_obj$log
  
}, error = function(e) {
  RSLT$mgwnbr_error <<- conditionMessage(e)
  cat(sprintf("[MGWNBR] Outer FAIL: %s\n",
              substr(conditionMessage(e), 1, 160)))
})

t_mg <- as.numeric((proc.time() - t0_mg)["elapsed"])

if (!is.null(m_mg)) {
  
  hp_mg <- m_mg$summary.hyperpar
  fe_mg <- m_mg$summary.fixed
  
  mu0_mg <- safe_fixed_mean(fe_mg, "mu0")
  mu1_mg <- safe_fixed_mean(fe_mg, "mu1")
  mu2_mg <- safe_fixed_mean(fe_mg, "mu2")
  
  rnd0_mg <- m_mg$summary.random$idx0
  rnd1_mg <- m_mg$summary.random$idx1
  rnd2_mg <- m_mg$summary.random$idx2
  
  b0_mg <- mu0_mg + rnd0_mg$mean
  b1_mg <- mu1_mg + rnd1_mg$mean
  b2_mg <- mu2_mg + rnd2_mg$mean
  
  b0_mg_lo <- mu0_mg + rnd0_mg$`0.025quant`
  b0_mg_hi <- mu0_mg + rnd0_mg$`0.975quant`
  
  b1_mg_lo <- mu1_mg + rnd1_mg$`0.025quant`
  b1_mg_hi <- mu1_mg + rnd1_mg$`0.975quant`
  
  b2_mg_lo <- mu2_mg + rnd2_mg$`0.025quant`
  b2_mg_hi <- mu2_mg + rnd2_mg$`0.975quant`
  
  err0_mg <- b0_mg - beta0_true
  err1_mg <- b1_mg - beta1_true
  err2_mg <- b2_mg - beta2_true
  
  RSLT$mgwnbr <- list(
    metrics = list(
      b0 = calc_metrics(b0_mg, beta0_true, b0_mg_lo, b0_mg_hi),
      b1 = calc_metrics(b1_mg, beta1_true, b1_mg_lo, b1_mg_hi),
      b2 = calc_metrics(b2_mg, beta2_true, b2_mg_lo, b2_mg_hi)
    ),
    
    bnd_rmse = list(
      b0_vh   = bnd_rmse(err0_mg, near_vhbnd),
      b1_diag = bnd_rmse(err1_mg, near_diag),
      b2_vh   = bnd_rmse(err2_mg, near_vhbnd),
      b2_diag = bnd_rmse(err2_mg, near_diag)
    ),
    
    bnd_profile = list(
      b0_vh   = bnd_profile(err0_mg, dist_to_vh),
      b1_diag = bnd_profile(err1_mg, dist_to_diag),
      b2_vh   = bnd_profile(err2_mg, dist_to_vh),
      b2_diag = bnd_profile(err2_mg, dist_to_diag)
    ),
    
    bandwidth = list(
      b0 = get_theta2(hp_mg, "idx0"),
      b1 = get_theta2(hp_mg, "idx1"),
      b2 = get_theta2(hp_mg, "idx2")
    ),
    
    ic = get_ic(m_mg),
    
    nb_size = get_nb_size(hp_mg),
    
    time = t_mg,
    
    fixed_ef = c(
      mu0 = mu0_mg,
      mu1 = mu1_mg,
      mu2 = mu2_mg
    ),
    
    fit_log = mg_fit_log
  )
  
  RSLT$fit_status$mgwnbr_ok <- TRUE
  
  cat(sprintf("[MGWNBR] OK | time=%.1fs | DIC=%.3f | WAIC=%.3f | LCPO=%.5f\n\n",
              t_mg,
              RSLT$mgwnbr$ic$dic,
              RSLT$mgwnbr$ic$waic,
              RSLT$mgwnbr$ic$lcpo))
  
  rm(m_mg, model.b0, model.b1, model.b2,
     rnd0_mg, rnd1_mg, rnd2_mg)
  gc(verbose = FALSE)
  
} else {
  
  if (is.null(RSLT$mgwnbr_error)) {
    RSLT$mgwnbr_error <- "all_mgwnbr_strategies_failed"
  }
  
  cat(sprintf("[MGWNBR] FAIL | time=%.1fs | error=%s\n\n",
              t_mg, RSLT$mgwnbr_error))
}

# ================================================================
# PART 9: MODEL 2 -- SVC-NBR with Matérn/SPDE fields
# ================================================================

t0_svc <- proc.time()

m_svc <- NULL

tryCatch({
  
  A_sim1 <- Diagonal(x = x1_rep) %*% A_sim_rep
  A_sim2 <- Diagonal(x = x2_rep) %*% A_sim_rep
  
  stk_svc <- inla.stack(
    data = list(y = y_rep),
    
    A = list(
      A_sim_rep,
      A_sim1,
      A_sim2,
      1
    ),
    
    effects = list(
      idx0_svc,
      idx1_svc,
      idx2_svc,
      list(
        mu0 = x0_rep,
        mu1 = x1_rep,
        mu2 = x2_rep
      )
    ),
    
    tag = "obs"
  )
  
  f_svc <- y ~ -1 + mu0 + mu1 + mu2 +
    f(
      s0,
      model = spde_sim_b0,
      extraconstr = mk_constr(spde_sim_b0)
    ) +
    f(
      s1,
      model = spde_sim_b1,
      extraconstr = mk_constr(spde_sim_b1)
    ) +
    f(
      s2,
      model = spde_sim_b2,
      extraconstr = mk_constr(spde_sim_b2)
    )
  
  nb_prior_svc <- list(
    hyper = list(
      theta = list(
        prior   = "normal",
        param   = c(log(NB_SIZE), 0.15),
        initial = log(NB_SIZE)
      )
    )
  )
  
  m_svc <- inla(
    f_svc,
    family = "nbinomial",
    data   = inla.stack.data(stk_svc),
    safe   = TRUE,
    
    control.family = nb_prior_svc,
    
    control.predictor = list(
      A       = inla.stack.A(stk_svc),
      compute = TRUE
    ),
    
    control.compute = list(
      dic  = TRUE,
      waic = TRUE,
      cpo  = TRUE,
      mlik = TRUE
    ),
    
    control.inla = list(
      strategy     = "simplified.laplace",
      int.strategy = "eb"
    ),
    
    verbose = FALSE
  )
  
}, error = function(e) {
  RSLT$svc_error <<- conditionMessage(e)
  cat(sprintf("[SVC-NBR] FAIL: %s\n",
              substr(conditionMessage(e), 1, 160)))
})

t_svc <- as.numeric((proc.time() - t0_svc)["elapsed"])

if (!is.null(m_svc)) {
  
  hp_svc <- m_svc$summary.hyperpar
  fe_svc <- m_svc$summary.fixed
  
  mu0_svc <- safe_fixed_mean(fe_svc, "mu0")
  mu1_svc <- safe_fixed_mean(fe_svc, "mu1")
  mu2_svc <- safe_fixed_mean(fe_svc, "mu2")
  
  rnd0_svc <- m_svc$summary.random$s0
  rnd1_svc <- m_svc$summary.random$s1
  rnd2_svc <- m_svc$summary.random$s2
  
  b0_svc <- mu0_svc + as.vector(A_sim %*% rnd0_svc$mean)
  b1_svc <- mu1_svc + as.vector(A_sim %*% rnd1_svc$mean)
  b2_svc <- mu2_svc + as.vector(A_sim %*% rnd2_svc$mean)
  
  # Approximate marginal interval projection.
  b0_svc_lo <- mu0_svc + as.vector(A_sim %*% rnd0_svc$`0.025quant`)
  b0_svc_hi <- mu0_svc + as.vector(A_sim %*% rnd0_svc$`0.975quant`)
  
  b1_svc_lo <- mu1_svc + as.vector(A_sim %*% rnd1_svc$`0.025quant`)
  b1_svc_hi <- mu1_svc + as.vector(A_sim %*% rnd1_svc$`0.975quant`)
  
  b2_svc_lo <- mu2_svc + as.vector(A_sim %*% rnd2_svc$`0.025quant`)
  b2_svc_hi <- mu2_svc + as.vector(A_sim %*% rnd2_svc$`0.975quant`)
  
  err0_svc <- b0_svc - beta0_true
  err1_svc <- b1_svc - beta1_true
  err2_svc <- b2_svc - beta2_true
  
  RSLT$svc <- list(
    metrics = list(
      b0 = calc_metrics(b0_svc, beta0_true, b0_svc_lo, b0_svc_hi),
      b1 = calc_metrics(b1_svc, beta1_true, b1_svc_lo, b1_svc_hi),
      b2 = calc_metrics(b2_svc, beta2_true, b2_svc_lo, b2_svc_hi)
    ),
    
    bnd_rmse = list(
      b0_vh   = bnd_rmse(err0_svc, near_vhbnd),
      b1_diag = bnd_rmse(err1_svc, near_diag),
      b2_vh   = bnd_rmse(err2_svc, near_vhbnd),
      b2_diag = bnd_rmse(err2_svc, near_diag)
    ),
    
    bnd_profile = list(
      b0_vh   = bnd_profile(err0_svc, dist_to_vh),
      b1_diag = bnd_profile(err1_svc, dist_to_diag),
      b2_vh   = bnd_profile(err2_svc, dist_to_vh),
      b2_diag = bnd_profile(err2_svc, dist_to_diag)
    ),
    
    hyperpar = list(
      b0 = extract_spde_hp(hp_svc, "s0"),
      b1 = extract_spde_hp(hp_svc, "s1"),
      b2 = extract_spde_hp(hp_svc, "s2")
    ),
    
    ic = get_ic(m_svc),
    
    nb_size = get_nb_size(hp_svc),
    
    time = t_svc,
    
    fixed_ef = c(
      mu0 = mu0_svc,
      mu1 = mu1_svc,
      mu2 = mu2_svc
    )
  )
  
  RSLT$fit_status$svc_ok <- TRUE
  
  cat(sprintf("[SVC-NBR] OK | time=%.1fs | DIC=%.3f | WAIC=%.3f | LCPO=%.5f\n\n",
              t_svc,
              RSLT$svc$ic$dic,
              RSLT$svc$ic$waic,
              RSLT$svc$ic$lcpo))
  
  rm(m_svc, rnd0_svc, rnd1_svc, rnd2_svc)
  gc(verbose = FALSE)
  
} else {
  
  if (is.null(RSLT$svc_error)) {
    RSLT$svc_error <- "svc_fit_failed"
  }
  
  cat(sprintf("[SVC-NBR] FAIL | time=%.1fs | error=%s\n\n",
              t_svc, RSLT$svc_error))
}

# ================================================================
# PART 10: Paired comparison
# ================================================================

RSLT$fit_status$both_ok <- isTRUE(RSLT$fit_status$mgwnbr_ok) &&
  isTRUE(RSLT$fit_status$svc_ok)

if (RSLT$fit_status$both_ok) {
  
  rmse_mg_b0  <- RSLT$mgwnbr$metrics$b0$rmse
  rmse_mg_b1  <- RSLT$mgwnbr$metrics$b1$rmse
  rmse_mg_b2  <- RSLT$mgwnbr$metrics$b2$rmse
  
  rmse_svc_b0 <- RSLT$svc$metrics$b0$rmse
  rmse_svc_b1 <- RSLT$svc$metrics$b1$rmse
  rmse_svc_b2 <- RSLT$svc$metrics$b2$rmse
  
  ratio_b0 <- rmse_mg_b0 / max(rmse_svc_b0, 1e-12)
  ratio_b1 <- rmse_mg_b1 / max(rmse_svc_b1, 1e-12)
  ratio_b2 <- rmse_mg_b2 / max(rmse_svc_b2, 1e-12)
  
  dic_mg   <- RSLT$mgwnbr$ic$dic
  waic_mg  <- RSLT$mgwnbr$ic$waic
  lcpo_mg  <- RSLT$mgwnbr$ic$lcpo
  
  dic_svc  <- RSLT$svc$ic$dic
  waic_svc <- RSLT$svc$ic$waic
  lcpo_svc <- RSLT$svc$ic$lcpo
  
  RSLT$comparison <- list(
    rmse_ratio = c(
      b0 = ratio_b0,
      b1 = ratio_b1,
      b2 = ratio_b2
    ),
    
    rmse_diff = c(
      b0 = rmse_mg_b0 - rmse_svc_b0,
      b1 = rmse_mg_b1 - rmse_svc_b1,
      b2 = rmse_mg_b2 - rmse_svc_b2
    ),
    
    rmse_win_mg = c(
      b0 = as.integer(rmse_mg_b0 < rmse_svc_b0),
      b1 = as.integer(rmse_mg_b1 < rmse_svc_b1),
      b2 = as.integer(rmse_mg_b2 < rmse_svc_b2)
    ),
    
    delta_ic = c(
      dic  = dic_mg  - dic_svc,
      waic = waic_mg - waic_svc,
      lcpo = lcpo_mg - lcpo_svc
    ),
    
    ic_win_mg = c(
      dic  = as.integer(is.finite(dic_mg)  && is.finite(dic_svc)  && dic_mg  < dic_svc),
      waic = as.integer(is.finite(waic_mg) && is.finite(waic_svc) && waic_mg < waic_svc),
      lcpo = as.integer(is.finite(lcpo_mg) && is.finite(lcpo_svc) && lcpo_mg < lcpo_svc)
    ),
    
    time_ratio = RSLT$mgwnbr$time / max(RSLT$svc$time, 1e-12)
  )
  
  cat("[Paired comparison]\n")
  cat(sprintf("  RMSE ratio MG/SVC: b0=%.3f | b1=%.3f | b2=%.3f\n",
              ratio_b0, ratio_b1, ratio_b2))
  cat(sprintf("  Delta IC MG-SVC : DIC=%.3f | WAIC=%.3f | LCPO=%.5f\n\n",
              RSLT$comparison$delta_ic["dic"],
              RSLT$comparison$delta_ic["waic"],
              RSLT$comparison$delta_ic["lcpo"]))
  
} else {
  
  RSLT$comparison <- NULL
}

save(RSLT, file = out_file)

# ================================================================
# ================================================================
#  Section 3: collect simulation results and generate plots 
#  collect simulation results from results_master
#  summary the results 
#  generate tables and plots 
# ================================================================
# ================================================================

setup_file  <- "sim_fixed_master.RData"
result_root <- "results_master"

RECORD_MISSING_IN_STATUS <- TRUE

if (!file.exists(setup_file)) {
  stop(sprintf("Cannot find setup file: %s", setup_file))
}

if (!dir.exists(result_root)) {
  stop(sprintf("Cannot find result directory: %s", result_root))
}

if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
}

load(setup_file)

scenario_ids <- scenario_table$scenario_id
R_expected   <- N_REP_MC

# ================================================================
# PART 1: load file
# ================================================================

load_RSLT_safe <- function(file_path) {
  
  e <- new.env(parent = emptyenv())
  
  ok <- tryCatch({
    load(file_path, envir = e)
    TRUE
  }, error = function(err) {
    warning(sprintf(
      "Failed to load %s | error: %s",
      file_path,
      conditionMessage(err)
    ))
    FALSE
  })
  
  if (!ok) {
    return(NULL)
  }
  
  if (!exists("RSLT", envir = e, inherits = FALSE)) {
    warning(sprintf("No object named RSLT in file: %s", file_path))
    return(NULL)
  }
  
  get("RSLT", envir = e, inherits = FALSE)
}

as_num <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(x[1]))
}

as_chr <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  as.character(x[1])
}

safe_status <- function(x, nm) {
  if (is.null(x) || !is.list(x)) return(FALSE)
  if (!nm %in% names(x)) return(FALSE)
  isTRUE(x[[nm]])
}

safe_value <- function(x, nm) {
  if (is.null(x)) return(NA_real_)
  
  if (is.list(x)) {
    if (!nm %in% names(x)) return(NA_real_)
    return(as_num(x[[nm]]))
  }
  
  if (!is.null(names(x)) && nm %in% names(x)) {
    return(as_num(x[[nm]]))
  }
  
  NA_real_
}

safe_metric <- function(x, nm) {
  safe_value(x, nm)
}

safe_ic <- function(x, nm) {
  safe_value(x, nm)
}

bind_rows_base <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}

q_fun <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

q025 <- function(x) q_fun(x, 0.025)
q25  <- function(x) q_fun(x, 0.25)
q50  <- function(x) q_fun(x, 0.50)
q75  <- function(x) q_fun(x, 0.75)
q975 <- function(x) q_fun(x, 0.975)

mcse_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_prop <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

get_model_obj <- function(RSLT, model_name) {
  if (is.null(RSLT) || !is.list(RSLT)) return(NULL)
  if (is.null(RSLT[[model_name]])) return(NULL)
  if (!is.list(RSLT[[model_name]])) return(NULL)
  RSLT[[model_name]]
}

normalise_RSLT <- function(RSLT, sid, kk, file_path) {
  
  if (is.null(RSLT) || !is.list(RSLT)) {
    return(list(
      fit_status = list(
        mgwnbr_ok = FALSE,
        svc_ok    = FALSE,
        both_ok   = FALSE
      ),
      mgwnbr_error = "RSLT_not_list",
      svc_error    = "RSLT_not_list"
    ))
  }
  
  if (is.null(RSLT$fit_status) || !is.list(RSLT$fit_status)) {
    RSLT$fit_status <- list(
      mgwnbr_ok = FALSE,
      svc_ok    = FALSE,
      both_ok   = FALSE
    )
  }
  
  if (!is.null(RSLT$mgwnbr) && !is.list(RSLT$mgwnbr)) {
    warning(sprintf(
      "RSLT$mgwnbr is atomic; set to NULL | scenario=%s | kk=%s | class=%s | file=%s",
      sid,
      kk,
      paste(class(RSLT$mgwnbr), collapse = ","),
      file_path
    ))
    
    RSLT$mgwnbr <- NULL
    RSLT$fit_status$mgwnbr_ok <- FALSE
    RSLT$mgwnbr_error <- "mgwnbr_atomic_or_malformed"
  }
  
  if (!is.null(RSLT$svc) && !is.list(RSLT$svc)) {
    warning(sprintf(
      "RSLT$svc is atomic; set to NULL | scenario=%s | kk=%s | class=%s | file=%s",
      sid,
      kk,
      paste(class(RSLT$svc), collapse = ","),
      file_path
    ))
    
    RSLT$svc <- NULL
    RSLT$fit_status$svc_ok <- FALSE
    RSLT$svc_error <- "svc_atomic_or_malformed"
  }
  
  RSLT$fit_status$mgwnbr_ok <- isTRUE(RSLT$fit_status$mgwnbr_ok) &&
    !is.null(get_model_obj(RSLT, "mgwnbr"))
  
  RSLT$fit_status$svc_ok <- isTRUE(RSLT$fit_status$svc_ok) &&
    !is.null(get_model_obj(RSLT, "svc"))
  
  RSLT$fit_status$both_ok <- isTRUE(RSLT$fit_status$mgwnbr_ok) &&
    isTRUE(RSLT$fit_status$svc_ok)
  
  RSLT
}

scenario_meta_row <- function(sid) {
  ss <- scenario_table[scenario_table$scenario_id == sid, , drop = FALSE]
  
  if (!nrow(ss)) {
    return(list(
      scenario_family = NA_character_,
      scenario_label  = NA_character_,
      delta           = NA_real_,
      delta_name      = NA_character_,
      SNR_b0          = NA_real_,
      SNR_b1          = NA_real_
    ))
  }
  
  list(
    scenario_family = as.character(ss$scenario_family[1]),
    scenario_label  = as.character(ss$scenario_label[1]),
    delta           = as.numeric(ss$delta[1]),
    delta_name      = as.character(ss$delta_name[1]),
    SNR_b0          = as.numeric(ss$SNR_b0[1]),
    SNR_b1          = as.numeric(ss$SNR_b1[1])
  )
}

make_meta_df <- function(sid, kk) {
  meta <- scenario_meta_row(sid)
  
  data.frame(
    scenario_id     = sid,
    kk              = kk,
    scenario_family = meta$scenario_family,
    scenario_label  = meta$scenario_label,
    delta           = meta$delta,
    delta_name      = meta$delta_name,
    SNR_b0          = meta$SNR_b0,
    SNR_b1          = meta$SNR_b1,
    stringsAsFactors = FALSE
  )
}

first_meta_from_df <- function(dd) {
  data.frame(
    scenario_id     = dd$scenario_id[1],
    scenario_family = dd$scenario_family[1],
    scenario_label  = dd$scenario_label[1],
    delta           = dd$delta[1],
    delta_name      = dd$delta_name[1],
    SNR_b0          = dd$SNR_b0[1],
    SNR_b1          = dd$SNR_b1[1],
    stringsAsFactors = FALSE
  )
}

split_apply <- function(df, group_vars, fun) {
  if (!is.data.frame(df) || !nrow(df)) return(data.frame())
  
  missing_vars <- setdiff(group_vars, names(df))
  if (length(missing_vars)) return(data.frame())
  
  key <- interaction(df[group_vars], drop = TRUE, lex.order = TRUE)
  groups <- split(df, key)
  out <- lapply(groups, fun)
  bind_rows_base(out)
}

# ================================================================
# PART 2: extract 
# ================================================================

extract_fit_status <- function(RSLT, sid, kk, file_path, file_exists = TRUE) {
  
  meta <- make_meta_df(sid, kk)
  mg_obj  <- get_model_obj(RSLT, "mgwnbr")
  svc_obj <- get_model_obj(RSLT, "svc")
  
  data.frame(
    meta,
    file_path    = file_path,
    file_exists  = file_exists,
    
    mgwnbr_ok    = safe_status(RSLT$fit_status, "mgwnbr_ok"),
    svc_ok       = safe_status(RSLT$fit_status, "svc_ok"),
    both_ok      = safe_status(RSLT$fit_status, "both_ok"),
    
    mgwnbr_error = ifelse(is.null(RSLT$mgwnbr_error),
                          NA_character_,
                          as_chr(RSLT$mgwnbr_error)),
    svc_error    = ifelse(is.null(RSLT$svc_error),
                          NA_character_,
                          as_chr(RSLT$svc_error)),
    
    mgwnbr_time  = if (!is.null(mg_obj)) as_num(mg_obj$time) else NA_real_,
    svc_time     = if (!is.null(svc_obj)) as_num(svc_obj$time) else NA_real_,
    
    stringsAsFactors = FALSE
  )
}

extract_global_metrics_one_model <- function(RSLT, sid, kk, model_name) {
  
  obj <- get_model_obj(RSLT, model_name)
  if (is.null(obj) || is.null(obj$metrics) || !is.list(obj$metrics)) {
    return(NULL)
  }
  
  coefs <- c("b0", "b1", "b2")
  out <- list()
  
  for (cc in coefs) {
    mm <- obj$metrics[[cc]]
    if (is.null(mm) || !is.list(mm)) next
    
    meta <- make_meta_df(sid, kk)
    
    out[[cc]] <- data.frame(
      meta,
      model       = model_name,
      coef        = cc,
      bias        = safe_metric(mm, "bias"),
      sd_err      = safe_metric(mm, "sd_err"),
      rmse        = safe_metric(mm, "rmse"),
      mae         = safe_metric(mm, "mae"),
      r2          = safe_metric(mm, "r2"),
      cov95       = safe_metric(mm, "cov95"),
      ci_width    = safe_metric(mm, "ci_width"),
      calib_ratio = safe_metric(mm, "calib_ratio"),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows_base(out)
}

extract_ic_one_model <- function(RSLT, sid, kk, model_name) {
  
  obj <- get_model_obj(RSLT, model_name)
  if (is.null(obj) || is.null(obj$ic) || !is.list(obj$ic)) {
    return(NULL)
  }
  
  meta <- make_meta_df(sid, kk)
  
  data.frame(
    meta,
    model = model_name,
    dic   = safe_ic(obj$ic, "dic"),
    waic  = safe_ic(obj$ic, "waic"),
    lcpo  = safe_ic(obj$ic, "lcpo"),
    mlik  = safe_ic(obj$ic, "mlik"),
    stringsAsFactors = FALSE
  )
}

extract_comparison_rmse <- function(RSLT, sid, kk) {
  
  if (is.null(RSLT$comparison) || !is.list(RSLT$comparison)) {
    return(NULL)
  }
  
  coefs <- c("b0", "b1", "b2")
  meta <- make_meta_df(sid, kk)
  
  out <- lapply(coefs, function(cc) {
    data.frame(
      meta,
      coef        = cc,
      rmse_ratio  = safe_value(RSLT$comparison$rmse_ratio, cc),
      rmse_diff   = safe_value(RSLT$comparison$rmse_diff, cc),
      rmse_win_mg = safe_value(RSLT$comparison$rmse_win_mg, cc),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows_base(out)
}

extract_comparison_ic <- function(RSLT, sid, kk) {
  
  if (is.null(RSLT$comparison) || !is.list(RSLT$comparison)) {
    return(NULL)
  }
  
  criteria <- c("dic", "waic", "lcpo")
  meta <- make_meta_df(sid, kk)
  
  out <- lapply(criteria, function(ic) {
    data.frame(
      meta,
      criterion = ic,
      delta_ic  = safe_value(RSLT$comparison$delta_ic, ic),
      ic_win_mg = safe_value(RSLT$comparison$ic_win_mg, ic),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows_base(out)
}

extract_boundary_profile_one_model <- function(RSLT, sid, kk, model_name) {
  
  obj <- get_model_obj(RSLT, model_name)
  if (is.null(obj) || is.null(obj$bnd_profile) || !is.list(obj$bnd_profile)) {
    return(NULL)
  }
  
  meta <- make_meta_df(sid, kk)
  profile_names <- names(obj$bnd_profile)
  out <- list()
  
  for (pn in profile_names) {
    
    prof <- obj$bnd_profile[[pn]]
    if (is.null(prof)) next
    
    parts <- strsplit(pn, "_")[[1]]
    coef_name <- parts[1]
    boundary_name <- ifelse(length(parts) >= 2L, parts[2], NA_character_)
    
    bins <- names(prof)
    if (is.null(bins) || length(bins) != length(prof)) {
      if (exists("BND_LABELS")) {
        bins <- BND_LABELS[seq_along(prof)]
      } else {
        bins <- paste0("bin", seq_along(prof))
      }
    }
    
    out[[pn]] <- data.frame(
      meta,
      model    = model_name,
      coef     = coef_name,
      boundary = boundary_name,
      bin      = bins,
      rmse     = as.numeric(prof),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows_base(out)
}

extract_count_diag <- function(RSLT, sid, kk) {
  
  if (is.null(RSLT$count_diag)) {
    return(NULL)
  }
  
  cd <- RSLT$count_diag
  meta <- make_meta_df(sid, kk)
  
  data.frame(
    meta,
    count_mean     = safe_value(cd, "mean"),
    count_var      = safe_value(cd, "var"),
    count_var_mean = safe_value(cd, "var_mean"),
    count_p_zero   = safe_value(cd, "p_zero"),
    count_q50      = safe_value(cd, "q50"),
    count_q90      = safe_value(cd, "q90"),
    count_q95      = safe_value(cd, "q95"),
    count_q99      = safe_value(cd, "q99"),
    count_max      = safe_value(cd, "max"),
    stringsAsFactors = FALSE
  )
}

extract_bandwidth_mgwnbr <- function(RSLT, sid, kk) {
  
  obj <- get_model_obj(RSLT, "mgwnbr")
  if (is.null(obj) || is.null(obj$bandwidth) || !is.list(obj$bandwidth)) {
    return(NULL)
  }
  
  meta <- make_meta_df(sid, kk)
  coefs <- names(obj$bandwidth)
  
  out <- lapply(coefs, function(cc) {
    
    bw <- obj$bandwidth[[cc]]
    if (is.null(bw) || !is.list(bw)) bw <- list()
    
    data.frame(
      meta,
      model   = "mgwnbr",
      coef    = cc,
      k_mean  = safe_value(bw, "k_mean"),
      k_q025  = safe_value(bw, "k_q025"),
      k_q975  = safe_value(bw, "k_q975"),
      bw_mean = safe_value(bw, "bw_mean"),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows_base(out)
}

extract_spde_hyperpar <- function(RSLT, sid, kk) {
  
  obj <- get_model_obj(RSLT, "svc")
  if (is.null(obj) || is.null(obj$hyperpar) || !is.list(obj$hyperpar)) {
    return(NULL)
  }
  
  meta <- make_meta_df(sid, kk)
  coefs <- names(obj$hyperpar)
  
  out <- lapply(coefs, function(cc) {
    
    hp <- obj$hyperpar[[cc]]
    if (is.null(hp) || !is.list(hp)) hp <- list()
    
    data.frame(
      meta,
      model       = "svc",
      coef        = cc,
      range_mean  = safe_value(hp, "range_mean"),
      range_q025  = safe_value(hp, "range_q025"),
      range_q975  = safe_value(hp, "range_q975"),
      sigma_mean  = safe_value(hp, "sigma_mean"),
      sigma_q025  = safe_value(hp, "sigma_q025"),
      sigma_q975  = safe_value(hp, "sigma_q975"),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows_base(out)
}

# ================================================================
# PART 3: Loop over
# ================================================================

raw_fit_status_list       <- list()
raw_global_metrics_list   <- list()
raw_ic_by_model_list      <- list()
raw_comparison_rmse_list  <- list()
raw_comparison_ic_list    <- list()
raw_boundary_profile_list <- list()
raw_count_diag_list       <- list()
raw_bandwidth_list        <- list()
raw_spde_hyperpar_list    <- list()

counter <- 1L

for (sid in scenario_ids) {
  
  cat(sprintf("[Collect] Scenario: %s\n", sid))
  
  scenario_dir <- file.path(result_root, sid)
  
  if (!dir.exists(scenario_dir)) {
    
    warning(sprintf("Scenario directory does not exist: %s", scenario_dir))
    
    if (RECORD_MISSING_IN_STATUS) {
      for (kk_miss in seq_len(R_expected)) {
        
        missing_file_path <- file.path(
          scenario_dir,
          sprintf("%d.RData", kk_miss)
        )
        
        empty_RSLT <- list(
          fit_status = list(
            mgwnbr_ok = FALSE,
            svc_ok    = FALSE,
            both_ok   = FALSE
          ),
          mgwnbr_error = "missing_scenario_dir",
          svc_error    = "missing_scenario_dir"
        )
        
        raw_fit_status_list[[counter]] <- extract_fit_status(
          empty_RSLT,
          sid,
          kk_miss,
          missing_file_path,
          file_exists = FALSE
        )
        
        counter <- counter + 1L
      }
    }
    
    next
  }
  
  files_this <- list.files(
    scenario_dir,
    pattern = "^[0-9]+\\.RData$",
    full.names = TRUE
  )
  
  if (length(files_this) == 0L) {
    
    warning(sprintf("No RData files found for scenario: %s", sid))
    
    if (RECORD_MISSING_IN_STATUS) {
      for (kk_miss in seq_len(R_expected)) {
        
        missing_file_path <- file.path(
          scenario_dir,
          sprintf("%d.RData", kk_miss)
        )
        
        empty_RSLT <- list(
          fit_status = list(
            mgwnbr_ok = FALSE,
            svc_ok    = FALSE,
            both_ok   = FALSE
          ),
          mgwnbr_error = "missing_file",
          svc_error    = "missing_file"
        )
        
        raw_fit_status_list[[counter]] <- extract_fit_status(
          empty_RSLT,
          sid,
          kk_miss,
          missing_file_path,
          file_exists = FALSE
        )
        
        counter <- counter + 1L
      }
    }
    
    next
  }
  
  kk_this <- suppressWarnings(
    as.integer(sub("\\.RData$", "", basename(files_this)))
  )
  
  keep <- is.finite(kk_this)
  files_this <- files_this[keep]
  kk_this <- kk_this[keep]
  
  ord <- order(kk_this)
  files_this <- files_this[ord]
  kk_this <- kk_this[ord]
  
  expected_kk <- seq_len(R_expected)
  missing_kk <- setdiff(expected_kk, kk_this)
  
  cat(sprintf(
    "  Existing result files: %d / %d\n",
    length(files_this),
    R_expected
  ))
  
  if (length(missing_kk) > 0L) {
    cat(sprintf(
      "  Missing files: %d. First missing kk: %s\n",
      length(missing_kk),
      paste(head(missing_kk, 20), collapse = ", ")
    ))
    
    if (length(missing_kk) > 20L) {
      cat("  ...\n")
    }
  }
  
  if (RECORD_MISSING_IN_STATUS && length(missing_kk) > 0L) {
    
    for (kk_miss in missing_kk) {
      
      missing_file_path <- file.path(
        scenario_dir,
        sprintf("%d.RData", kk_miss)
      )
      
      empty_RSLT <- list(
        fit_status = list(
          mgwnbr_ok = FALSE,
          svc_ok    = FALSE,
          both_ok   = FALSE
        ),
        mgwnbr_error = "missing_file",
        svc_error    = "missing_file"
      )
      
      raw_fit_status_list[[counter]] <- extract_fit_status(
        empty_RSLT,
        sid,
        kk_miss,
        missing_file_path,
        file_exists = FALSE
      )
      
      counter <- counter + 1L
    }
  }
  
  for (ii in seq_along(files_this)) {
    
    kk <- kk_this[ii]
    file_path <- files_this[ii]
    
    cat(sprintf("  Loading existing file: kk=%d | %s\n",
                kk, basename(file_path)))
    
    RSLT <- load_RSLT_safe(file_path)
    
    if (is.null(RSLT)) {
      
      bad_RSLT <- list(
        fit_status = list(
          mgwnbr_ok = FALSE,
          svc_ok    = FALSE,
          both_ok   = FALSE
        ),
        mgwnbr_error = "load_failed_or_missing_RSLT",
        svc_error    = "load_failed_or_missing_RSLT"
      )
      
      raw_fit_status_list[[counter]] <- extract_fit_status(
        bad_RSLT,
        sid,
        kk,
        file_path,
        file_exists = TRUE
      )
      
      counter <- counter + 1L
      next
    }
    
    RSLT <- normalise_RSLT(RSLT, sid, kk, file_path)
    
    raw_fit_status_list[[counter]] <- extract_fit_status(
      RSLT,
      sid,
      kk,
      file_path,
      file_exists = TRUE
    )
    
    raw_global_metrics_list[[counter]] <- bind_rows_base(list(
      extract_global_metrics_one_model(RSLT, sid, kk, "mgwnbr"),
      extract_global_metrics_one_model(RSLT, sid, kk, "svc")
    ))
    
    raw_ic_by_model_list[[counter]] <- bind_rows_base(list(
      extract_ic_one_model(RSLT, sid, kk, "mgwnbr"),
      extract_ic_one_model(RSLT, sid, kk, "svc")
    ))
    
    raw_comparison_rmse_list[[counter]] <- extract_comparison_rmse(
      RSLT,
      sid,
      kk
    )
    
    raw_comparison_ic_list[[counter]] <- extract_comparison_ic(
      RSLT,
      sid,
      kk
    )
    
    raw_boundary_profile_list[[counter]] <- bind_rows_base(list(
      extract_boundary_profile_one_model(RSLT, sid, kk, "mgwnbr"),
      extract_boundary_profile_one_model(RSLT, sid, kk, "svc")
    ))
    
    raw_count_diag_list[[counter]] <- extract_count_diag(
      RSLT,
      sid,
      kk
    )
    
    raw_bandwidth_list[[counter]] <- extract_bandwidth_mgwnbr(
      RSLT,
      sid,
      kk
    )
    
    raw_spde_hyperpar_list[[counter]] <- extract_spde_hyperpar(
      RSLT,
      sid,
      kk
    )
    
    counter <- counter + 1L
  }
}

# ================================================================
# PART 4: Bind raw objects
# ================================================================

raw_fit_status        <- bind_rows_base(raw_fit_status_list)
raw_global_metrics    <- bind_rows_base(raw_global_metrics_list)
raw_ic_by_model       <- bind_rows_base(raw_ic_by_model_list)
raw_comparison_rmse   <- bind_rows_base(raw_comparison_rmse_list)
raw_comparison_ic     <- bind_rows_base(raw_comparison_ic_list)
raw_boundary_profile  <- bind_rows_base(raw_boundary_profile_list)
raw_count_diag        <- bind_rows_base(raw_count_diag_list)
raw_bandwidth_mgwnbr  <- bind_rows_base(raw_bandwidth_list)
raw_spde_hyperpar     <- bind_rows_base(raw_spde_hyperpar_list)

# ================================================================
# PART 5: Summary, tables and figures
# ================================================================
################ Figure 1 #######################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(sf)
library(grid)
library(scales)

plot_df_list <- lapply(names(beta_true_list), function(sid) {
  bt <- beta_true_list[[sid]]
  data.frame(
    scenario_id = sid,
    scenario_family = bt$scenario_family,
    scenario_label = bt$scenario_label,
    delta = bt$delta,
    delta_name = bt$delta_name,
    unit_id = seq_len(N_LOC),
    east = coords[, 1],
    north = coords[, 2],
    beta0 = bt$beta0,
    beta1 = bt$beta1,
    beta2 = bt$beta2
  )
})

plot_df <- bind_rows(plot_df_list)

# Make ordered scenario labels
plot_df$scenario_id <- factor(
  plot_df$scenario_id,
  levels = scenario_table$scenario_id
)

plot_long <- plot_df %>%
  pivot_longer(
    cols = c(beta0, beta1, beta2),
    names_to = "coef",
    values_to = "value"
  ) %>%
  mutate(
    coef = factor(coef, levels = c("beta0", "beta1", "beta2")),
    scenario_id = factor(scenario_id, levels = scenario_table$scenario_id)
  )

p_true <- ggplot(plot_long, aes(x = east, y = north, fill = value)) +
  geom_raster() +
  facet_grid(coef ~ scenario_id) +
  coord_equal() +
  scale_fill_viridis_c(option = "viridis") +
  labs(
    #title = "True coefficient surfaces under different simulation settings",
    x = "East",
    y = "North",
    fill = "Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 9, face = "bold"),
    strip.text.y = element_text(size = 9, face = "bold")
  )

p_true

################ Figure 2 #######################################
summary_rmse_ratio <- split_apply(
  raw_comparison_rmse,
  group_vars = c("scenario_id", "coef"),
  fun = function(dd) {
    
    meta <- first_meta_from_df(dd)
    
    rr  <- dd$rmse_ratio
    rd  <- dd$rmse_diff
    win <- dd$rmse_win_mg
    
    data.frame(
      meta,
      coef = dd$coef[1],
      
      n_both = sum(is.finite(rr)),
      
      mean_ratio   = safe_mean(rr),
      sd_ratio     = safe_sd(rr),
      mcse_ratio   = mcse_mean(rr),
      median_ratio = safe_median(rr),
      q025_ratio   = q025(rr),
      q25_ratio    = q25(rr),
      q75_ratio    = q75(rr),
      q975_ratio   = q975(rr),
      
      mean_rmse_diff   = safe_mean(rd),
      median_rmse_diff = safe_median(rd),
      
      mg_rmse_win_prob = safe_prop(win),
      
      winner_by_mean_ratio = ifelse(is.finite(safe_mean(rr)) &&
                                      safe_mean(rr) < 1,
                                    "MGWNBR", "SVC-NBR"),
      winner_by_median_ratio = ifelse(is.finite(safe_median(rr)) &&
                                        safe_median(rr) < 1,
                                      "MGWNBR", "SVC-NBR"),
      
      stringsAsFactors = FALSE
    )
  }
)

if (nrow(summary_rmse_ratio)) {
  
  summary_rmse_ratio$scenario_id <- factor(
    summary_rmse_ratio$scenario_id,
    levels = scenario_levels
  )
  
  summary_rmse_ratio$scenario_label_plot <- scenario_labels[
    as.character(summary_rmse_ratio$scenario_id)
  ]
  
  summary_rmse_ratio$scenario_label_plot <- factor(
    summary_rmse_ratio$scenario_label_plot,
    levels = scenario_labels[scenario_levels]
  )
  
  summary_rmse_ratio$coef <- factor(
    summary_rmse_ratio$coef,
    levels = coef_levels
  )
  
  summary_rmse_ratio$winner <- ifelse(
    summary_rmse_ratio$mean_ratio < 1,
    "MGWNBR",
    "SVC-NBR"
  )
  
  summary_rmse_ratio$log2_mean_ratio <- log2(summary_rmse_ratio$mean_ratio)
  summary_rmse_ratio$log2_median_ratio <- log2(summary_rmse_ratio$median_ratio)
  
  summary_rmse_ratio$ratio_label <- fmt_num(summary_rmse_ratio$mean_ratio, 2)
}

p1 <- ggplot(
  summary_rmse_ratio,
  aes(
    x = scenario_label_plot,
    y = coef,
    fill = log2_mean_ratio
  )
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = ratio_label),
    size = 3.5,
    color = "black"
  ) +
  scale_y_discrete(
    labels = coef_labels
  ) +
  scale_fill_gradient2(
    name = "log(RMSE ratio)\nMGWNBR / SVC-NBR",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -max(abs(summary_rmse_ratio$log2_mean_ratio), na.rm = TRUE),
      max(abs(summary_rmse_ratio$log2_mean_ratio), na.rm = TRUE)
    ),
    oob = scales::squish
  ) +
  labs(
    title = "Global coefficient recovery: RMSE ratio",
    subtitle = "Blue indicates lower RMSE for MGWNBR; red indicates lower RMSE for SVC-NBR",
    x = "Scenario",
    y = "Coefficient"
  ) +
  theme_pub(12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p1


######### Table 3 ###############################
summary_ic <- split_apply(
  raw_comparison_ic,
  group_vars = c("scenario_id", "criterion"),
  fun = function(dd) {
    
    meta <- first_meta_from_df(dd)
    
    dx  <- dd$delta_ic
    win <- dd$ic_win_mg
    
    data.frame(
      meta,
      criterion = dd$criterion[1],
      
      n_both = sum(is.finite(dx)),
      
      mean_delta_ic   = safe_mean(dx),
      sd_delta_ic     = safe_sd(dx),
      mcse_delta_ic   = mcse_mean(dx),
      median_delta_ic = safe_median(dx),
      q025_delta_ic   = q025(dx),
      q25_delta_ic    = q25(dx),
      q75_delta_ic    = q75(dx),
      q975_delta_ic   = q975(dx),
      
      mg_ic_win_prob = safe_prop(win),
      
      winner_by_median_delta = ifelse(is.finite(safe_median(dx)) &&
                                        safe_median(dx) < 0,
                                      "MGWNBR", "SVC-NBR"),
      
      stringsAsFactors = FALSE
    )
  }
)

############### Figure 3 ##################################
if (nrow(summary_ic)) {
  
  summary_ic$scenario_id <- factor(
    summary_ic$scenario_id,
    levels = scenario_levels
  )
  
  summary_ic$scenario_label_plot <- scenario_labels[
    as.character(summary_ic$scenario_id)
  ]
  
  summary_ic$scenario_label_plot <- factor(
    summary_ic$scenario_label_plot,
    levels = scenario_labels[scenario_levels]
  )
  
  summary_ic$criterion <- factor(
    summary_ic$criterion,
    levels = criterion_levels
  )
  
  summary_ic$criterion_label <- criterion_labels[
    as.character(summary_ic$criterion)
  ]
  
  summary_ic$winner <- ifelse(
    summary_ic$median_delta_ic < 0,
    "MGWNBR",
    "SVC-NBR"
  )
  
  # For heatmap color, use signed log transform to reduce extreme IC values.
  # Negative means MGWNBR is better.
  summary_ic$signed_log_delta <- sign(summary_ic$median_delta_ic) *
    log10(abs(summary_ic$median_delta_ic) + 1)
  
  summary_ic$delta_label <- ifelse(
    abs(summary_ic$median_delta_ic) >= 100,
    formatC(summary_ic$median_delta_ic, format = "f", digits = 0),
    formatC(summary_ic$median_delta_ic, format = "f", digits = 2)
  )
}

p2 <- ggplot(
  summary_ic,
  aes(
    x = scenario_label_plot,
    y = criterion,
    fill = signed_log_delta
  )
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = delta_label),
    size = 3.3,
    color = "black"
  ) +
  scale_y_discrete(
    labels = criterion_labels
  ) +
  scale_fill_gradient2(
    name = expression(paste("signed", log(Delta(IC)))),
    #name = "signed log10\n|median Delta IC|",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    title = "Information-criterion comparison",
    subtitle = expression(Delta * "IC = IC"[MGWNBR] - "IC"[SVC] * "; blue favours MGWNBR, red favours SVC-NBR"),
    x = "Scenario",
    y = "Criterion"
  ) +
  theme_pub(12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p2

################################## Figure 4 ##########################################
summary_boundary_ratio <- split_apply(
  raw_boundary_ratio,
  group_vars = c("scenario_id", "coef", "boundary", "bin"),
  fun = function(dd) {
    
    meta <- first_meta_from_df(dd)
    
    rr  <- dd$rmse_ratio
    win <- dd$mg_win
    
    data.frame(
      meta,
      coef     = dd$coef[1],
      boundary = dd$boundary[1],
      bin      = dd$bin[1],
      
      n_both = sum(is.finite(rr)),
      
      mean_ratio   = safe_mean(rr),
      sd_ratio     = safe_sd(rr),
      mcse_ratio   = mcse_mean(rr),
      median_ratio = safe_median(rr),
      q025_ratio   = q025(rr),
      q25_ratio    = q25(rr),
      q75_ratio    = q75(rr),
      q975_ratio   = q975(rr),
      
      mg_win_prob = safe_prop(win),
      
      winner_by_mean_ratio = ifelse(is.finite(safe_mean(rr)) &&
                                      safe_mean(rr) < 1,
                                    "MGWNBR", "SVC-NBR"),
      winner_by_median_ratio = ifelse(is.finite(safe_median(rr)) &&
                                        safe_median(rr) < 1,
                                      "MGWNBR", "SVC-NBR"),
      
      stringsAsFactors = FALSE
    )
  }
)

if (nrow(summary_boundary_ratio)) {
  
  summary_boundary_ratio$scenario_id <- factor(
    summary_boundary_ratio$scenario_id,
    levels = scenario_levels
  )
  
  summary_boundary_ratio$scenario_label_plot <- scenario_labels[
    as.character(summary_boundary_ratio$scenario_id)
  ]
  
  summary_boundary_ratio$scenario_label_plot <- factor(
    summary_boundary_ratio$scenario_label_plot,
    levels = scenario_labels[scenario_levels]
  )
  
  summary_boundary_ratio$coef <- factor(
    summary_boundary_ratio$coef,
    levels = coef_levels
  )
  
  summary_boundary_ratio$bin <- factor(
    summary_boundary_ratio$bin,
    levels = bin_levels
  )
  
  summary_boundary_ratio$boundary_label <- boundary_labels[
    summary_boundary_ratio$boundary
  ]
  
  summary_boundary_ratio$boundary_label <- ifelse(
    is.na(summary_boundary_ratio$boundary_label),
    summary_boundary_ratio$boundary,
    summary_boundary_ratio$boundary_label
  )
  
  summary_boundary_ratio$coef_boundary <- paste0(
    as.character(summary_boundary_ratio$coef),
    " / ",
    summary_boundary_ratio$boundary
  )
  
  summary_boundary_ratio$coef_boundary <- factor(
    summary_boundary_ratio$coef_boundary,
    levels = c("b0 / vh", "b1 / diag", "b2 / vh", "b2 / diag")
  )
}

bd_main <- summary_boundary_ratio[
  summary_boundary_ratio$scenario_family == "H" &
    summary_boundary_ratio$coef_boundary %in% c("b0 / vh", "b1 / diag"),
]

bd_main$coef_boundary_label <- ifelse(
  bd_main$coef_boundary == "b0 / vh",
  "beta0 near vertical/horizontal boundary",
  "beta1 near diagonal boundary"
)

bd_main$coef_boundary_label <- factor(
  bd_main$coef_boundary_label,
  levels = c(
    "beta0 near vertical/horizontal boundary",
    "beta1 near diagonal boundary"
  )
)

p3 <- ggplot(
  bd_main,
  aes(
    x = bin,
    y = median_ratio,
    group = scenario_label_plot,
    color = scenario_label_plot
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey30"
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~ coef_boundary_label,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 2.5),
    labels = c("0.5", "0.75", "1", "1.25", "1.5", "2", "2.5")
  ) +
  scale_color_manual(
    name = "Hard scenario",
    values = c(
      "Hard d=0.05" = "#A6CEE3",
      "Hard d=0.10" = "#1F78B4",
      "Hard d=0.15" = "#B2DF8A",
      "Hard d=0.25" = "#33A02C"
    )
  ) +
  labs(
    title = "Boundary-distance RMSE ratio under hard discontinuities",
    subtitle = "MGWNBR improves recovery near discontinuity boundaries, while SVC-NBR remains better away from boundaries",
    x = "Distance to discontinuity boundary",
    y = "Median RMSE ratio: MGWNBR / SVC-NBR"
  ) +
  theme_pub(12) +
  theme(
    legend.position = "top"
  )

p3

############################## Figure 5 ##########################################
bd_win <- summary_boundary_ratio[
  summary_boundary_ratio$scenario_family == "H" &
    summary_boundary_ratio$coef_boundary %in% c("b0 / vh", "b1 / diag"),
]

bd_win$coef_boundary_label <- ifelse(
  bd_win$coef_boundary == "b0 / vh",
  "beta0 near vertical/horizontal boundary",
  "beta1 near diagonal boundary"
)

bd_win$coef_boundary_label <- factor(
  bd_win$coef_boundary_label,
  levels = c(
    "beta0 near vertical/horizontal boundary",
    "beta1 near diagonal boundary"
  )
)

p4 <- ggplot(
  bd_win,
  aes(
    x = bin,
    y = mg_win_prob,
    group = scenario_label_plot,
    color = scenario_label_plot
  )
) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey30"
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~ coef_boundary_label,
    ncol = 1
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    name = "Hard scenario",
    values = c(
      "Hard d=0.05" = "#A6CEE3",
      "Hard d=0.10" = "#1F78B4",
      "Hard d=0.15" = "#B2DF8A",
      "Hard d=0.25" = "#33A02C"
    )
  ) +
  labs(
    title = "Boundary-distance MGWNBR win probability",
    subtitle = "Win probability is highest near hard discontinuity boundaries",
    x = "Distance to discontinuity boundary",
    y = "Pr(MGWNBR RMSE < SVC-NBR RMSE)"
  ) +
  theme_pub(12) +
  theme(
    legend.position = "top"
  )

p4