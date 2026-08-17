# ==============================================================
# Chicago Crime 2023: Empirical Analysis
# Comparing MGWNBR vs SVC-NBR via INLA
# ==============================================================
# Data: chicago_boundary_count_model_2023.csv
# Y:     crime_count
# x1:    log population density
# x2:    violent crime proportion
# Offset: log_population_offset
# Boundary: Police District hard boundary (same_side_matrix.csv)
# ==============================================================

library(INLA); library(Matrix); library(MASS)
library(sf); library(dplyr); library(FNN)
library(ggplot2); library(viridis); library(patchwork)
library(spdep)

inla.setOption(smtp= "taucs")
inla.setOption(num.threads = "4:1")

set.seed(2025)

# ================================================================
# SECTION 1: load data
# ================================================================

cat("=== Loading Chicago Crime Data ===\n")

dat_raw <- read.csv(
  "chicago_boundary_count_model_2023.csv",
  stringsAsFactors = FALSE
)
same_side_raw <- read.csv(
  "same_side_matrix.csv",
  row.names = 1, check.names = FALSE
)

dat <- dat_raw %>%
  filter(model_include_recommended == 1) %>%
  arrange(analysis_unit_id)

N_LOC <- nrow(dat)            
cat(sprintf("N_LOC = %d\n", N_LOC))

uid<- as.character(dat$analysis_unit_id)
same_side_mat <- as.matrix(same_side_raw[uid, uid])
diag(same_side_mat) <- 1L

# ── covariates ─────────────────────────────────────────────────
dat <- dat %>% mutate(
  pop_density = population_est_area_weighted / pmax(community_area_km2, 1e-4),
  x1_raw= log(pop_density + 1),            
  x2_raw      = violent_crime_count / pmax(crime_count, 1), 
  x1_std      = as.numeric(scale(x1_raw)),
  x2_std      = as.numeric(scale(x2_raw))
)

x0_vec<- rep(1.0, N_LOC)    
mu0_vec <- rep(1.0, N_LOC)    
mu1_vec <- dat$x1_std          
idx_loc <- seq_len(N_LOC)     

# ── overdispersion ───────────────────────────────────────────────
mu_y <- mean(dat$crime_count)
vr_y <- var(dat$crime_count)
cat(sprintf("Overdispersion: mean=%.1f  var=%.1f  V/M=%.2f\n",
            mu_y, vr_y, vr_y / mu_y))

# NB_SIZE
NB_SIZE <- max(mu_y^2 / (vr_y - mu_y), 0.5)
NB_SIZE <- min(NB_SIZE, 100.0)
cat(sprintf("NB_SIZE prior centre: %.2f\n", NB_SIZE))

# ================================================================
# SECTION 2: projection & kNN
# ================================================================

cat("\n=== Building Spatial Objects ===\n")

# projection
coords_sf<- st_as_sf(dat, coords = c("longitude","latitude"), crs = 4326)
coords_proj <- st_coordinates(st_transform(coords_sf, 32616))
coords_km   <- coords_proj / 1000.0
dat$x_km    <- coords_km[, 1]
dat$y_km    <- coords_km[, 2]

# kNN
K_PRE <- min(40L, N_LOC - 1L)
K_MAX <- K_PRE - 2L

knn_res<- FNN::get.knn(coords_km, k = K_PRE)
knn_idx_mat   <- knn_res$nn.index   # N_LOC × K_PRE
knn_dist_raw  <- knn_res$nn.dist    # N_LOC × K_PRE

d_max_g<- max(knn_dist_raw)
D_knn_mat     <- knn_dist_raw / d_max_g   

med_bw_at_k_vec <- apply(D_knn_mat, 2, median) * d_max_g  

cat(sprintf("K_PRE=%d  K_MAX=%d  d_max=%.2f km\n", K_PRE, K_MAX, d_max_g))

# hyper-parameters
K0_CTR  <- 15L;K1_CTR  <- 18L;  K2_CTR  <- 12L
K_MIN_B0 <- 3L; K_MIN_B1 <- 3L; K_MIN_B2 <- 3L
B0_KSHAPE <- 4.0; B1_KSHAPE <- 4.0; B2_KSHAPE <- 5.0
K0_RATE   <- K0_CTR * (B0_KSHAPE +1)
K1_RATE   <- K1_CTR * (B1_KSHAPE + 1)
K2_RATE   <- K2_CTR * (B2_KSHAPE + 1)

V_SIGMA_B0 <- 1.5; V_SIGMA_B1 <- 0.6; V_SIGMA_B2 <- 0.6
LAMBDA_B0  <- -log(0.5) / sqrt(V_SIGMA_B0)
LAMBDA_B1  <- -log(0.5) / sqrt(V_SIGMA_B1)
LAMBDA_B2  <- -log(0.5) / sqrt(V_SIGMA_B2)
TH1_B0     <- log(1.0 / V_SIGMA_B0)
TH1_B1     <- log(1.0 / V_SIGMA_B1)
TH1_B2     <- log(1.0 / V_SIGMA_B2)

# ================================================================
# SECTION 3: boundary indicator
# ================================================================
near_bnd <- apply(same_side_mat == 0L, 1, any)
cat(sprintf("Near-boundary units: %d / %d (%.1f%%)\n",
            sum(near_bnd), N_LOC, 100 * mean(near_bnd)))

dist_to_bnd <- sapply(seq_len(N_LOC), function(i) {
  cross_j <- which(same_side_mat[i, ] == 0L)
  if (!length(cross_j)) return(Inf)
  min(sqrt((dat$x_km[i] - dat$x_km[cross_j])^2 +(dat$y_km[i] - dat$y_km[cross_j])^2)) / 2.0
})
dist_to_bnd[!is.finite(dist_to_bnd)]<-
  max(dist_to_bnd[is.finite(dist_to_bnd)])

# ================================================================
# SECTION 4: SPDE 
# ================================================================
mesh_hb <- inla.mesh.2d(
  loc      = coords_km,
  max.edge = c(5.0, 15.0),
  offset   = c(4.0, 10.0),
  cutoff   = 2.0
)
cat(sprintf("SPDE mesh: %d nodes\n", mesh_hb$n))

A_hb <- inla.spde.make.A(mesh_hb, loc = coords_km)   # N_LOC × mesh_n

spde_hb_b0 <- inla.spde2.pcmatern(mesh_hb, alpha = 2,prior.range = c(8.0, 0.05), prior.sigma = c(sqrt(V_SIGMA_B0), 0.05))
spde_hb_b1 <- inla.spde2.pcmatern(mesh_hb, alpha = 2,
                                  prior.range = c(8.0, 0.05), prior.sigma = c(sqrt(V_SIGMA_B1), 0.05))
spde_hb_b2 <- inla.spde2.pcmatern(mesh_hb, alpha = 2,
                                  prior.range = c(8.0, 0.05), prior.sigma = c(sqrt(V_SIGMA_B2), 0.05))

idx0_hb <- inla.spde.make.index("s0", n.spde = mesh_hb$n)
idx1_hb <- inla.spde.make.index("s1", n.spde = mesh_hb$n)
idx2_hb <- inla.spde.make.index("s2", n.spde = mesh_hb$n)

# ================================================================
# SECTION 5: rgeneric
# ================================================================

inla.rgeneric.mgwnbr.harddgp <- function(cmd= c("graph","Q","mu","initial","log.norm.const","log.prior","quit"),
                                         theta = NULL) {
  envir <- parent.env(environment())
  itp <- function() {
    tau<- exp(max(min(theta[1L], 12.0), -4.0))
    k_raw <- exp(max(min(theta[2L], log(K_MAX)), log(K_MIN)))
    list(tau = tau, k = k_raw)
  }
  
  graph <- function() {
    K_G<- min(K_PRE_G, 20L)
    i_idx <- rep(seq_len(n_g), each = K_G)
    j_idx <- as.vector(t(KNN_IDX_G[, seq_len(K_G)]))
    all_i <- c(i_idx, j_idx, seq_len(n_g))
    all_j <- c(j_idx, i_idx, seq_len(n_g))
    
    G <- Matrix::sparseMatrix(
      i            = all_i,
      j            = all_j,
      x            = rep(1.0, length(all_i)),
      dims         = c(n_g, n_g),
      use.last.ij  = TRUE
    )
    G
  }
  
  Q <- function() {
    p<- itp()
    k_lo <- max(1L, floor(p$k))
    k_hi <- min(K_PRE_G, ceiling(p$k))
    al   <- p$k - k_lo
    
    bw_lo <- KNN_D_G[k_lo, ]
    bw_hi <- KNN_D_G[k_hi, ]
    b_vec <- (1.0 - al) * bw_lo + al * bw_hi
    b_vec <- pmax(b_vec, 1e-4)
    
    K_USE <- min(
      K_PRE_G,
      max(ceiling(p$k) + 5L, 8L),
      floor(n_g *0.15)
    )
    K_USE <- max(K_USE, 5L)
    
    i_idx  <- rep(seq_len(n_g), each = K_USE)
    j_idx  <- as.vector(t(KNN_IDX_G[, seq_len(K_USE)]))
    d_vals <- as.vector(t(KNN_D_G[seq_len(K_USE), ]))
    b_src  <- rep(b_vec, each = K_USE)
    w_vals <- exp(-d_vals^2 / (2.0 * b_src^2))
    w_vals[!is.finite(w_vals)] <- 0.0
    w_vals[w_vals < 1e-7]<- 0.0
    
    W_sp <- Matrix::sparseMatrix(
      i    = i_idx,
      j    = j_idx,
      x    = w_vals,
      dims = c(n_g, n_g)
    )
    W_sym <- (W_sp + Matrix::t(W_sp)) / 2.0
    rs    <- Matrix::rowSums(W_sym)
    nug   <- p$tau * max(rs, 1.0) * 5e-2+ 1e-5
    
    Q_out <- p$tau * (Matrix::Diagonal(x = rs) - W_sym) +nug   *  Matrix::Diagonal(n_g)
    
    methods::as(Q_out, "dgCMatrix")
  }
  
  mu<- function() rep(0.0, n_g)
  log.norm.const <- function() numeric(0)
  
  log.prior <- function() {
    tau <- exp(theta[1L])
    if (!is.finite(tau) || tau <= 0) return(-1e10)
    sigma<- 1.0 / sqrt(tau)
    lp_tau <- log(LAMBDA_G /2) - theta[1L] / 2 - LAMBDA_G * sigma
    lp_k<- -(K_SHAPE_G + 1) * theta[2L] - K_RATE_G * exp(-theta[2L])
    lp     <- lp_tau + lp_k
    if (!is.finite(lp)) return(-1e10)
    lp
  }
  
  initial <- function() c(TH1_INIT_G, log(K_CTR_G))
  quit<- function() invisible()
  
  if (!length(theta) || is.null(theta)) theta <- initial()
  do.call(match.arg(cmd), args = list())
}

# ================================================================
# SECTION 6: additional functions
# ================================================================

mk_constr   <- function(sp) list(A = matrix(1, 1, sp$n.spde), e = 0)
mk_constr_n <- function(nn) list(A = matrix(1, 1, nn), e = 0)

bnd_spread <- function(coef_vec, mask)list(near= sd(coef_vec[ mask]),
                                           far   = sd(coef_vec[!mask]),
                                           ratio = sd(coef_vec[mask]) / max(sd(coef_vec[!mask]), 1e-8))

bnd_profile <- function(coef_vec, dist_vec,
                        breaks = c(0, 1, 2, 4, 8, Inf)) {
  cuts <- cut(dist_vec, breaks = breaks, right = FALSE,
              labels = paste0("d<", breaks[-1]))
  tapply(coef_vec, cuts, sd, na.rm = TRUE)
}

extract_spde_hp <- function(hp, fname) {
  rn<- rownames(hp)
  r_row <- grep(sprintf("Range.*%s", fname), rn, ignore.case = TRUE)
  s_row <- grep(sprintf("Stdev.*%s", fname), rn, ignore.case = TRUE)
  out   <- list()
  if (length(r_row)) {
    out$range_mean <- hp[r_row[1], "mean"]
    out$range_q025 <- hp[r_row[1], "0.025quant"]
    out$range_q975 <- hp[r_row[1], "0.975quant"]
  }
  if (length(s_row)) out$sigma_mean <- hp[s_row[1], "mean"]
  out
}

get_theta2 <- function(hp, tag) {
  row <- grep(sprintf("Theta2.*%s", tag), rownames(hp), ignore.case = TRUE)
  if (!length(row))
    return(list(k_mean=NA_real_, k_q025=NA_real_,
                k_q975=NA_real_, bw_mean_km=NA_real_))
  k_m <- exp(hp[row[1], "mean"])
  ki<- max(1L, min(round(k_m), length(med_bw_at_k_vec)))
  list(k_mean     = k_m,
       k_q025     = exp(hp[row[1], "0.025quant"]),
       k_q975     = exp(hp[row[1], "0.975quant"]),
       bw_mean_km = med_bw_at_k_vec[ki])
}

get_lcpo <- function(fit) {
  cpo_v <- fit$cpo$cpo
  valid <- cpo_v > 0 & is.finite(cpo_v) & !is.na(cpo_v)
  if (!sum(valid)) NA_real_ else -mean(log(cpo_v[valid]))
}

mk_model_hb <- function(th1, k_ctr, k_min, k_shape, k_rate, lam) {
  inla.rgeneric.define(
    inla.rgeneric.mgwnbr.harddgp,
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

fit_mgwnbr_robust <- function(f, data_df, nb_prior, verbose = TRUE) {
  strategies <- list(
    list(tag= "Gaussian+EB",
         ci   = list(strategy="gaussian", int.strategy="eb",
                     diagonal=1.0, h=0.01, tolerance=1e-3),
         full = FALSE),
    list(tag  = "SimpLaplace+EB",
         ci   = list(strategy="simplified.laplace", int.strategy="eb",
                     diagonal=1e-2, h=8e-3, tolerance=1e-4),
         full = TRUE),
    list(tag  = "SimpLaplace+EB(safe)",
         ci   = list(strategy="simplified.laplace", int.strategy="eb",
                     diagonal=0.1, h=0.02, tolerance=1e-3),
         full = TRUE)
  )
  theta_mode <- NULL; m_final <- NULL
  for (s in strategies) {
    cat(sprintf("  -> [%s]...\n", s$tag))
    ctrl_mode <- if (!is.null(theta_mode))
      list(theta=theta_mode, restart=TRUE) else list(result=NULL)
    m_try <- tryCatch(
      inla(f, family="nbinomial", data=data_df, safe=TRUE,
           control.family= nb_prior,
           control.predictor = list(compute=s$full),
           control.compute   = if (s$full)
             list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE) else
               list(dic=FALSE, waic=FALSE, cpo=FALSE, mlik=FALSE),
           control.mode= ctrl_mode,
           control.inla  = s$ci,
           verbose       = verbose),
      error = function(e) {
        cat(sprintf("     FAIL: %s\n", substr(conditionMessage(e), 1, 80)))
        NULL
      }
    )
    if (!is.null(m_try)) {
      theta_mode <- m_try$mode$theta
      cat(sprintf("     OK | theta=%s\n",
                  paste(round(theta_mode, 2), collapse =",")))
      if (s$full) { m_final <- m_try; break }
    }
  }
  m_final
}

# ================================================================
# SECTION 7: MODEL1-Bayesian MGWNBR
# ================================================================

model.b0 <- mk_model_hb(TH1_B0, K0_CTR, K_MIN_B0, B0_KSHAPE, K0_RATE, LAMBDA_B0)
model.b1 <- mk_model_hb(TH1_B1, K1_CTR, K_MIN_B1, B1_KSHAPE, K1_RATE, LAMBDA_B1)
model.b2 <- mk_model_hb(TH1_B2, K2_CTR, K_MIN_B2, B2_KSHAPE, K2_RATE, LAMBDA_B2)

# data_mg
data_mg <- data.frame(
  y= dat$crime_count,
  x0= x0_vec,        
  x1      = dat$x1_std,
  x2      = dat$x2_std,
  mu0     = mu0_vec,
  mu1     = mu1_vec,         
  idx0    = idx_loc,
  idx1    = idx_loc,
  idx2    = idx_loc,
  log_pop = dat$log_population_offset
)

f_mg <- y ~ -1 + mu0 + mu1 +
  f(idx0, x0, model=model.b0, n=N_LOC, extraconstr=mk_constr_n(N_LOC)) +
  f(idx1, x1, model=model.b1, n=N_LOC, extraconstr=mk_constr_n(N_LOC)) +
  f(idx2, x2, model=model.b2, n=N_LOC) +
  offset(log_pop)                

nb_prior_mg <- list(hyper = list(theta = list(
  prior= "normal",
  param= c(log(NB_SIZE),0.20),
  initial = log(NB_SIZE)
)))

t0_mg <- proc.time()
m_mg  <- fit_mgwnbr_robust(f_mg, data_mg, nb_prior_mg)
t_mg  <- (proc.time() - t0_mg)["elapsed"]

RSLT_MG <- NULL
if (!is.null(m_mg)) {
  hp_mg  <- m_mg$summary.hyperpar
  fe_mg  <- m_mg$summary.fixed
  mu0_mg <- fe_mg["mu0", "mean"]
  mu1_mg <- fe_mg["mu1", "mean"]
  
  rnd0_mg <- m_mg$summary.random$idx0
  rnd1_mg <- m_mg$summary.random$idx1
  rnd2_mg <- m_mg$summary.random$idx2
  
  b0_mg    <- mu0_mg + rnd0_mg$mean
  b1_mg    <- mu1_mg + rnd1_mg$mean
  b2_mg    <- rnd2_mg$mean
  b0_mg_lo <- mu0_mg + rnd0_mg$`0.025quant`
  b0_mg_hi <- mu0_mg + rnd0_mg$`0.975quant`
  b1_mg_lo <- mu1_mg + rnd1_mg$`0.025quant`
  b1_mg_hi <- mu1_mg + rnd1_mg$`0.975quant`
  b2_mg_lo <- rnd2_mg$`0.025quant`
  b2_mg_hi <- rnd2_mg$`0.975quant`
  
  nb_row_mg <- grep("size|nbinom", rownames(hp_mg), ignore.case=TRUE)[1]
  
  RSLT_MG <- list(
    coef = data.frame(
      b0=b0_mg, b0_lo=b0_mg_lo, b0_hi=b0_mg_hi,
      b1=b1_mg, b1_lo=b1_mg_lo, b1_hi=b1_mg_hi,
      b2=b2_mg, b2_lo=b2_mg_lo, b2_hi=b2_mg_hi
    ),
    bnd_spread = list(
      b0 = bnd_spread(b0_mg, near_bnd),
      b1 = bnd_spread(b1_mg, near_bnd)
    ),
    bnd_profile = list(
      b0 = bnd_profile(b0_mg, dist_to_bnd),
      b1 = bnd_profile(b1_mg, dist_to_bnd)
    ),
    bandwidth = list(
      bw_b0 = get_theta2(hp_mg, "idx0"),
      bw_b1 = get_theta2(hp_mg, "idx1"),
      bw_b2 = get_theta2(hp_mg, "idx2")
    ),
    ic = list(
      dic= m_mg$dic$dic,
      waic = m_mg$waic$waic,
      lcpo = get_lcpo(m_mg),
      mlik = as.numeric(m_mg$mlik[1, 1])
    ),
    nb_size= hp_mg[nb_row_mg, "mean"],
    time     = t_mg,
    fixed_ef = c(mu0=mu0_mg, mu1=mu1_mg)
  )
  
  cat(sprintf("[MGWNBR] OK | %.1fs | DIC=%.2f | LCPO=%.4f | NB_size=%.2f\n",
              t_mg, m_mg$dic$dic, RSLT_MG$ic$lcpo, RSLT_MG$nb_size))
  rm(m_mg, model.b0, model.b1, model.b2, rnd0_mg, rnd1_mg, rnd2_mg)
  gc(verbose = TRUE)
} else {
  cat(sprintf("[MGWNBR] FAIL | %.1fs\n", t_mg))
}

# ================================================================
# SECTION 8: MODEL 2-SVC-NBR
# ================================================================

A_hb1 <- Diagonal(x = dat$x1_std) %*% A_hb
A_hb2 <- Diagonal(x = dat$x2_std) %*% A_hb

stk_hb <- inla.stack(
  data    = list(y = dat$crime_count),
  A       = list(A_hb, A_hb1, A_hb2,1),
  effects = list(
    idx0_hb, idx1_hb, idx2_hb,
    list(mu0     = mu0_vec,
         mu1     = dat$x1_std,
         mu2     = dat$x2_std,
         log_pop = dat$log_population_offset)
  ),
  tag = "obs"
)

f_svc_hb <- y ~ -1 + mu0 + mu1 + mu2 +
  f(s0, model=spde_hb_b0, extraconstr=mk_constr(spde_hb_b0)) +
  f(s1, model=spde_hb_b1, extraconstr=mk_constr(spde_hb_b1)) +
  f(s2, model=spde_hb_b2, extraconstr=mk_constr(spde_hb_b2)) +
  offset(log_pop)                

nb_prior_svc <- list(hyper = list(theta = list(
  prior   = "normal",
  param   = c(log(NB_SIZE), 0.20),
  initial = log(NB_SIZE)
)))

t0_svc <- proc.time()
m_svc  <- tryCatch(
  inla(f_svc_hb,family          = "nbinomial",
       data            = inla.stack.data(stk_hb),
       safe            = TRUE,
       control.family  = nb_prior_svc,
       control.predictor = list(A=inla.stack.A(stk_hb), compute=TRUE),
       control.compute   = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE),
       control.inla      = list(strategy="simplified.laplace",int.strategy="eb"),
       verbose = TRUE),
  error = function(e) {
    cat(sprintf("[SVC-NBR] FAIL: %s\n", substr(conditionMessage(e),1,80)))
    NULL
  }
)
t_svc <- (proc.time() - t0_svc)["elapsed"]

RSLT_SVC <- NULL
if (!is.null(m_svc)) {
  hp_svc  <- m_svc$summary.hyperpar
  fe_svc  <- m_svc$summary.fixed
  mu0_svc <- fe_svc["mu0", "mean"]
  mu1_svc <- fe_svc["mu1", "mean"]
  mu2_svc <- fe_svc["mu2", "mean"]
  
  rnd0_svc <- m_svc$summary.random$s0
  rnd1_svc <- m_svc$summary.random$s1
  rnd2_svc <- m_svc$summary.random$s2
  
  b0_svc    <- mu0_svc + as.vector(A_hb %*% rnd0_svc$mean)
  b1_svc    <- mu1_svc + as.vector(A_hb %*% rnd1_svc$mean)
  b2_svc    <- mu2_svc + as.vector(A_hb %*% rnd2_svc$mean)
  b0_svc_lo <- mu0_svc + as.vector(A_hb %*% rnd0_svc$`0.025quant`)
  b0_svc_hi <- mu0_svc + as.vector(A_hb %*% rnd0_svc$`0.975quant`)
  b1_svc_lo <- mu1_svc + as.vector(A_hb %*% rnd1_svc$`0.025quant`)
  b1_svc_hi <- mu1_svc + as.vector(A_hb %*% rnd1_svc$`0.975quant`)
  b2_svc_lo <- mu2_svc + as.vector(A_hb %*% rnd2_svc$`0.025quant`)
  b2_svc_hi <- mu2_svc + as.vector(A_hb %*% rnd2_svc$`0.975quant`)
  
  nb_row_svc <- grep("size|nbinom", rownames(hp_svc), ignore.case=TRUE)[1]
  
  RSLT_SVC <- list(
    coef = data.frame(
      b0=b0_svc, b0_lo=b0_svc_lo, b0_hi=b0_svc_hi,
      b1=b1_svc, b1_lo=b1_svc_lo, b1_hi=b1_svc_hi,
      b2=b2_svc, b2_lo=b2_svc_lo, b2_hi=b2_svc_hi
    ),
    bnd_spread = list(
      b0 = bnd_spread(b0_svc, near_bnd),
      b1 = bnd_spread(b1_svc, near_bnd)
    ),
    bnd_profile = list(
      b0 = bnd_profile(b0_svc, dist_to_bnd),
      b1 = bnd_profile(b1_svc, dist_to_bnd)
    ),
    hyperpar = list(
      hp_b0 = extract_spde_hp(hp_svc, "s0"),
      hp_b1 = extract_spde_hp(hp_svc, "s1"),
      hp_b2 = extract_spde_hp(hp_svc, "s2")
    ),
    ic = list(
      dic  = m_svc$dic$dic,
      waic = m_svc$waic$waic,
      lcpo = get_lcpo(m_svc),
      mlik = as.numeric(m_svc$mlik[1, 1])
    ),
    nb_size  = hp_svc[nb_row_svc, "mean"],
    time     = t_svc,
    fixed_ef = c(mu0=mu0_svc, mu1=mu1_svc, mu2=mu2_svc)
  )
  
  cat(sprintf("[SVC-NBR] OK | %.1fs | DIC=%.2f | LCPO=%.4f\n",
              t_svc, m_svc$dic$dic, RSLT_SVC$ic$lcpo))
}

# ================================================================
# SECTION 9: outputs
# ================================================================

#### Table 4: model comparison #####

ic_tab <- data.frame(
  Model= c("MGWNBR", "SVC-NBR"),
  DIC     = c(RSLT_MG$ic$dic,RSLT_SVC$ic$dic),
  WAIC    = c(RSLT_MG$ic$waic,RSLT_SVC$ic$waic),
  negLCPO = c(RSLT_MG$ic$lcpo,  RSLT_SVC$ic$lcpo),
  NB_size = c(RSLT_MG$nb_size,  RSLT_SVC$nb_size),
  Time_s  = c(RSLT_MG$time,     RSLT_SVC$time)
)
ic_tab$dDIC  <- round(ic_tab$DIC  - min(ic_tab$DIC),2)
ic_tab$dWAIC <- round(ic_tab$WAIC - min(ic_tab$WAIC), 2)
print(ic_tab)

res_tab <- data.frame(
  unit_id      = dat$analysis_unit_id,
  longitude    = dat$longitude,
  latitude     = dat$latitude,
  crime_count  = dat$crime_count,
  near_bnd     = near_bnd,
  dist_bnd_km  = dist_to_bnd,
  b0_mg=RSLT_MG$coef$b0, b0_mg_lo=RSLT_MG$coef$b0_lo, b0_mg_hi=RSLT_MG$coef$b0_hi,
  b1_mg=RSLT_MG$coef$b1, b2_mg=RSLT_MG$coef$b2,
  b0_svc=RSLT_SVC$coef$b0, b1_svc=RSLT_SVC$coef$b1, b2_svc=RSLT_SVC$coef$b2,
  b0_diff = abs(RSLT_MG$coef$b0 - RSLT_SVC$coef$b0)
)


dist_q<- quantile(dist_to_bnd, probs = c(0.25, 0.5, 0.75))
bnd_zone <- cut(dist_to_bnd,
                breaks = c(0, dist_q[1], dist_q[2], dist_q[3], Inf),
                labels = c("Q1", "Q2", "Q3", "Q4"),
                include.lowest = TRUE)

spread_by_zone <- function(coef_vec) {
  tapply(coef_vec, bnd_zone, sd, na.rm = TRUE)
}


# ================================================================
# Boundary-distance discrepancy profile
# ================================================================

cat("\n=== Boundary-distance discrepancy profile ===\n")

if (exists("dist_to_cpd_bnd_km")) {
  dist_boundary_km <- dist_to_cpd_bnd_km
  dist_source <- "sf_exact_cpd_boundary"
} else {
  dist_boundary_km <- dist_to_bnd
  dist_source <- "approx_cross_district_half_distance"
}

cat(sprintf("[Distance source] %s\n", dist_source))

# ------------------------------------------------
# 1. Compute absolute inter-model discrepancies
# ------------------------------------------------

delta_b0 <- abs(RSLT_MG$coef$b0 - RSLT_SVC$coef$b0)
delta_b1 <- abs(RSLT_MG$coef$b1 - RSLT_SVC$coef$b1)
delta_b2 <- abs(RSLT_MG$coef$b2 - RSLT_SVC$coef$b2)

# Optional signed differences, useful for diagnostics
signed_delta_b0 <- RSLT_MG$coef$b0 - RSLT_SVC$coef$b0
signed_delta_b1 <- RSLT_MG$coef$b1 - RSLT_SVC$coef$b1
signed_delta_b2 <- RSLT_MG$coef$b2 - RSLT_SVC$coef$b2

# ------------------------------------------------
# 2. Distance bins
# ------------------------------------------------

DIST_BREAKS_KM <- c(0, 0.5, 1, 2, 4, 8, Inf)
DIST_LABELS_KM <- c("[0,0.5)", "[0.5,1)", "[1,2)", "[2,4)", "[4,8)", "[8,Inf)")

dist_bin_fixed <- cut(
  dist_boundary_km,
  breaks = DIST_BREAKS_KM,
  labels = DIST_LABELS_KM,
  right = FALSE,
  include.lowest = TRUE
)

dist_q <- quantile(
  dist_boundary_km,
  probs = c(0, 0.25, 0.50, 0.75, 1),
  na.rm = TRUE
)

dist_q_unique <- unique(as.numeric(dist_q))

if (length(dist_q_unique) >= 3) {
  dist_bin_quantile <- cut(
    dist_boundary_km,
    breaks = dist_q_unique,
    include.lowest = TRUE,
    right = FALSE
  )
} else {
  dist_bin_quantile <- dist_bin_fixed
}

USE_BIN <- "fixed" 

if (USE_BIN == "fixed") {
  dist_bin_use <- dist_bin_fixed
  bin_type <- "Fixed km bins"
} else {
  dist_bin_use <- dist_bin_quantile
  bin_type <- "Distance quantile bins"
}

# ------------------------------------------------
# 3. Unit-level discrepancy data
# ------------------------------------------------

boundary_discrepancy_df <- data.frame(
  unit_id = dat$analysis_unit_id,
  longitude = dat$longitude,
  latitude = dat$latitude,
  crime_count = dat$crime_count,
  log_crime_count = log1p(dat$crime_count),
  
  dist_boundary_km = dist_boundary_km,
  dist_bin = dist_bin_use,
  dist_bin_fixed = dist_bin_fixed,
  dist_bin_quantile = dist_bin_quantile,
  near_bnd = near_bnd,
  
  b0_mg = RSLT_MG$coef$b0,
  b1_mg = RSLT_MG$coef$b1,
  b2_mg = RSLT_MG$coef$b2,
  
  b0_svc = RSLT_SVC$coef$b0,
  b1_svc = RSLT_SVC$coef$b1,
  b2_svc = RSLT_SVC$coef$b2,
  
  delta_b0 = delta_b0,
  delta_b1 = delta_b1,
  delta_b2 = delta_b2,
  
  signed_delta_b0 = signed_delta_b0,
  signed_delta_b1 = signed_delta_b1,
  signed_delta_b2 = signed_delta_b2
)

# ------------------------------------------------
# 4. Summary function by distance bin
# ------------------------------------------------

se_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

summ_one_delta <- function(df, delta_col) {
  df %>%
    group_by(dist_bin) %>%
    summarise(
      n = n(),
      mean_delta = mean(.data[[delta_col]], na.rm = TRUE),
      median_delta = median(.data[[delta_col]], na.rm = TRUE),
      sd_delta = sd(.data[[delta_col]], na.rm = TRUE),
      se_delta = se_mean(.data[[delta_col]]),
      q25_delta = quantile(.data[[delta_col]], 0.25, na.rm = TRUE),
      q75_delta = quantile(.data[[delta_col]], 0.75, na.rm = TRUE),
      mean_dist_km = mean(dist_boundary_km, na.rm = TRUE),
      median_dist_km = median(dist_boundary_km, na.rm = TRUE),
      .groups = "drop"
    )
}

profile_b0 <- summ_one_delta(boundary_discrepancy_df, "delta_b0") %>%
  mutate(coef = "b0")

profile_b1 <- summ_one_delta(boundary_discrepancy_df, "delta_b1") %>%
  mutate(coef = "b1")

profile_b2 <- summ_one_delta(boundary_discrepancy_df, "delta_b2") %>%
  mutate(coef = "b2")

boundary_discrepancy_profile <- bind_rows(profile_b0, profile_b1, profile_b2) %>%
  mutate(
    coef = factor(coef, levels = c("b0", "b1", "b2")),
    dist_bin = factor(dist_bin, levels = levels(dist_bin_use))
  ) %>%
  arrange(coef, dist_bin)


# ------------------------------------------------
# 5. Spearman correlation tests
# ------------------------------------------------

cor_b0 <- cor.test(
  boundary_discrepancy_df$dist_boundary_km,
  boundary_discrepancy_df$delta_b0,
  method = "spearman",
  exact = FALSE
)

cor_b1 <- cor.test(
  boundary_discrepancy_df$dist_boundary_km,
  boundary_discrepancy_df$delta_b1,
  method = "spearman",
  exact = FALSE
)

cor_b2 <- cor.test(
  boundary_discrepancy_df$dist_boundary_km,
  boundary_discrepancy_df$delta_b2,
  method = "spearman",
  exact = FALSE
)

cor_tab <- data.frame(
  coef = c("b0", "b1", "b2"),
  spearman_rho = c(
    unname(cor_b0$estimate),
    unname(cor_b1$estimate),
    unname(cor_b2$estimate)
  ),
  p_value = c(
    cor_b0$p.value,
    cor_b1$p.value,
    cor_b2$p.value
  )
)

cat("\n[Spearman correlation: |Delta beta_k| vs distance to CPD boundary]\n")
print(cor_tab)

################### plots ##############################

library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(scales)

gdf_community <- st_read("chicago_community_areas_current.geojson",
                         quiet = TRUE)
gdf_districts <- st_read("chicago_police_districts_current.geojson",
                         quiet = TRUE)

gdf_community <- st_transform(gdf_community, 4326)
gdf_districts <- st_transform(gdf_districts, 4326)

cat("\n=== Computing exact distance to CPD district boundaries ===\n")

gdf_districts_utm <- st_transform(gdf_districts, 32616)

cpd_boundary_lines_utm <- gdf_districts_utm |>
  st_boundary() |>
  st_union()

unit_centroids_sf <- st_as_sf(
  dat,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

unit_centroids_utm <- st_transform(unit_centroids_sf, 32616)

dist_to_cpd_bnd_km <- as.numeric(st_distance(unit_centroids_utm, cpd_boundary_lines_utm)) / 1000

cat(sprintf("Distance to CPD boundary: min=%.3f km | median=%.3f km | max=%.3f km\n",
            min(dist_to_cpd_bnd_km, na.rm = TRUE),
            median(dist_to_cpd_bnd_km, na.rm = TRUE),
            max(dist_to_cpd_bnd_km, na.rm = TRUE)))


gdf_dist_boundary <- st_cast(
  st_union(gdf_districts) |> st_boundary(),
  "LINESTRING"
)
gdf_dist_lines <- gdf_districts |>
  st_cast("MULTILINESTRING") |>
  st_cast("LINESTRING")


coef_df <- data.frame(
  area_id= dat$community_area_num,        
  
  # MGWNBR
  b0_mg     = RSLT_MG$coef$b0,
  b1_mg     = RSLT_MG$coef$b1,
  b2_mg     = RSLT_MG$coef$b2,
  b0_mg_lo= RSLT_MG$coef$b0_lo,
  b0_mg_hi  = RSLT_MG$coef$b0_hi,
  b0_mg_ci  = RSLT_MG$coef$b0_hi - RSLT_MG$coef$b0_lo,
  
  # SVC-NBR
  b0_svc    = RSLT_SVC$coef$b0,
  b1_svc    = RSLT_SVC$coef$b1,
  b2_svc    = RSLT_SVC$coef$b2,
  b0_svc_ci = RSLT_SVC$coef$b0_hi - RSLT_SVC$coef$b0_lo,
  
  b0_diff   = abs(RSLT_MG$coef$b0 - RSLT_SVC$coef$b0),
  y_count= dat$crime_count
)

GEO_KEY <- "area_numbe"   

gdf_community[[GEO_KEY]] <- as.character(gdf_community[[GEO_KEY]])
coef_df$area_id           <- as.character(coef_df$area_id)

gdf_plot <- gdf_community |>
  left_join(coef_df, by = setNames("area_id", GEO_KEY))


# ── Theme ──────────────────────────────────────────────
theme_map <- theme_void(base_size = 11) +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 9,face = "bold"),
    legend.text        = element_text(size = 8),
    legend.key.height  = unit(1.2, "cm"),
    legend.key.width   = unit(0.4, "cm"),
    plot.title         = element_text(size = 11, face = "bold",
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 9, color = "grey40",
                                      margin = margin(b = 6)),
    plot.margin        = margin(6, 6, 6, 6)
  )

mk_choropleth <- function(gdf,
                          fill_col, 
                          title,      
                          subtitle= "",pal       = "plasma",    
                          limits = limits,
                          direction = 1,
                          legend_title = "Estimates",
                          dist_lines= gdf_dist_lines) {
  
  v_range <- range(gdf[[fill_col]], na.rm = TRUE)
  ggplot() +
    geom_sf(data   = gdf,
            aes(fill = .data[[fill_col]]),
            color  = "white",
            linewidth = 0.15) +
    geom_sf(data  = dist_lines,
            color = "black",
            linewidth = 0.7,
            alpha = 0.85,
            inherit.aes = FALSE) +
    scale_fill_viridis_c(
      option= pal,
      direction = direction,
      #limits    = c(-v_lim, v_lim),
      limits = limits,
      oob       = scales::squish,
      name      = legend_title,
      na.value  = "grey85",
      guide     = guide_colorbar(title.position = "top",
                                 barwidth = 0.8,
                                 barheight = 10)
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_map
}

##### Figure 8 (a)(d) ##########################
p_b0_mg  <- mk_choropleth(
  gdf_plot, "b0_mg",
  title=expression(paste("MGWNBR ", beta[0])),
  subtitle = sprintf("bandwidth ≈ %.1f km",
                     RSLT_MG$bandwidth$bw_b0$bw_mean_km),
  pal      = "plasma",
  limits = c(-4,-1)
)

p_b0_svc <- mk_choropleth(
  gdf_plot, "b0_svc",
  title=expression(paste("SVC-NBR ", beta[0])),
  subtitle = sprintf("spatial range ≈ %.1f km",
                     RSLT_SVC$hyperpar$hp_b0$range_mean),
  pal      = "plasma",
  limits = c(-4,-1)
)

##### Figure 8 (b)(e) ##########################
p_b1_mg  <- mk_choropleth(
  gdf_plot, "b1_mg",
  title    = expression(paste("MGWNBR ", beta[1])),
  subtitle = sprintf("bandwidth ≈ %.1f km",
                     RSLT_MG$bandwidth$bw_b1$bw_mean_km),
  pal      = "viridis",
  limits = c(-0.5,0.5)
)

p_b1_svc <- mk_choropleth(
  gdf_plot, "b1_svc",
  title    = expression(paste("SVC-NBR ", beta[1])),
  subtitle = sprintf("spatial range ≈ %.1f km",
                     RSLT_SVC$hyperpar$hp_b1$range_mean),
  pal      = "viridis",
  limits = c(-0.5,0.5)
)

##### Figure 8 (c)(f) ##########################
p_b2_mg  <- mk_choropleth(
  gdf_plot, "b2_mg",
  title    = expression(paste("MGWNBR ", beta[2])),
  subtitle = sprintf("bandwidth ≈ %.1f km",
                     RSLT_MG$bandwidth$bw_b2$bw_mean_km),
  pal      = "mako",
  limits = c(-2.5,2.5)
)

p_b2_svc <- mk_choropleth(
  gdf_plot, "b2_svc",
  title    = expression(paste("SVC-NBR ", beta[2])),
  subtitle = sprintf("spatial range ≈ %.1f km",
                     RSLT_SVC$hyperpar$hp_b2$range_mean),
  pal      = "mako",
  limits = c(-2.5,2.5)
)

############### Figure 9 ##################
p_diff <-ggplot() +
  geom_sf(data  = gdf_plot,
          aes(fill = b0_diff),
          color = "white", linewidth = 0.15) +
  geom_sf(data  = gdf_dist_lines,
          color = "black", linewidth = 0.9, alpha = 0.9) +
  scale_fill_viridis_c(
    option = "plasma", name = "value",
    na.value = "grey85",
    guide = guide_colorbar(title.position = "top",barwidth = 0.8, barheight = 10)
  ) +
  labs(title= expression(paste(abs(Delta(beta[0])), "  (MGWNBR − SVC-NBR)")))+
  theme_map

############# Figure 7 #########################
p_crime <- ggplot() +
  geom_sf(data  = gdf_plot,
          aes(fill = log1p(y_count)),
          color = "white", linewidth = 0.15) +
  geom_sf(data  = gdf_dist_lines,
          color = "black", linewidth = 0.9, alpha = 0.9) +
  scale_fill_viridis_c(
    option   = "rocket",
    #name     = "log(crime counts+1)",
    name = "value",
    na.value = "grey85",
    guide    = guide_colorbar(title.position = "top",
                              barwidth = 0.8, barheight = 10)
  ) +
  labs(title    = "crime counts (log scale)") +
  theme_map


################# Figure 11 #########################

profile_all_plot <- boundary_discrepancy_profile %>%
  mutate(
    coef_label = recode(
      as.character(coef),
      b0 = "Intercept: beta[0]",
      b1 = "Population density: beta[1]",
      b2 = "Violent-crime composition: beta[2]"
    ),
    ymin = pmax(mean_delta - 1.96 * se_delta, 0),
    ymax = mean_delta + 1.96 * se_delta,
    dist_bin = factor(dist_bin, levels = levels(dist_bin_use))
  )

p_delta_all_profile <- ggplot(
  profile_all_plot,
  aes(x = dist_bin, y = mean_delta, group = coef, color = coef, fill = coef)
) +
  geom_ribbon(
    aes(ymin = ymin, ymax = ymax),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  facet_wrap(
    ~coef_label,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(
      b0 = "#08519C",
      b1 = "#238B45",
      b2 = "#B2182B"
    ),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(
      b0 = "#08519C",
      b1 = "#238B45",
      b2 = "#B2182B"
    ),
    guide = "none"
  ) +
  labs(
    title = expression(paste("Boundary-distance profiles of ", "|", Delta, beta[k], "|")),
    subtitle = "Absolute MGWNBR--SVC-NBR coefficient discrepancies by distance to CPD district boundary",
    x = "Distance to nearest CPD district boundary (km)",
    y = expression(paste("Mean ", "|", Delta, beta[k], "|"))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey70"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey35"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  )

print(p_delta_all_profile)


#################### Figure 10 ##################

p_delta_b0_scatter <- ggplot(
  boundary_discrepancy_df,
  aes(x = dist_boundary_km, y = delta_b0)
) +
  geom_point(
    aes(color = near_bnd),
    alpha = 0.70,
    size = 2.3
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    color = "#B2182B",
    fill = "#F4A582",
    linewidth = 1.2,
    alpha = 0.20
  ) +
  scale_color_manual(
    values = c("TRUE" = "#D95F0E", "FALSE" = "#3182BD"),
    labels = c("TRUE" = "Near cross-district unit", "FALSE" = "Other unit")
  ) +
  labs(
    title = expression(paste("|", Delta, beta[0], "| vs. distance to CPD boundary")),
    subtitle = sprintf(
      "Spearman rho = %.3f, p = %.3g",
      cor_tab$spearman_rho[cor_tab$coef == "b0"],
      cor_tab$p_value[cor_tab$coef == "b0"]
    ),
    x = "Distance to nearest CPD district boundary (km)",
    y = expression(paste("|", Delta, beta[0], "|")),
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey70"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey35"),
    legend.position = "top"
  )

print(p_delta_b0_scatter)


