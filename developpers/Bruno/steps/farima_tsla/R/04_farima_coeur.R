# =============================================================================
# 04_farima_coeur.R  --  Noyau mathématique du FARIMA(p, d, q)
# =============================================================================
# Briques de bas niveau, SANS estimation, toutes testées dans tests/test_modules.R :
#   - coefficients du différenciateur fractionnaire (1 - B)^d et de son inverse
#   - poids MA(infini) et AR(infini)
#   - densité spectrale
#   - autocovariances EXACTES (convolution FI(d) x ARMA)
#   - algorithme de Durbin-Levinson (innovations, vraisemblance exacte)
#   - simulation EXACTE (Davies-Harte)
#   - reparamétrisation garantissant stationnarité (|d|<0.5, racines AR) et inversibilité (MA)
# =============================================================================

# ---- 1. Séries binomiales de (1-B)^d et (1-B)^(-d) ----------------------------------------------
# (1 - B)^d = sum_j pi_j B^j,  pi_0 = 1,  pi_j = pi_{j-1} (j - 1 - d) / j
coef_diff_frac <- function(d, n) {
  if (n <= 1) return(1)
  j <- seq_len(n - 1)
  cumprod(c(1, (j - 1 - d) / j))
}
# (1 - B)^(-d) = sum_j psi_j B^j,  psi_0 = 1,  psi_j = psi_{j-1} (j - 1 + d) / j   ~  j^(d-1)/Gamma(d)
coef_ma_frac <- function(d, n) {
  if (n <= 1) return(1)
  j <- seq_len(n - 1)
  cumprod(c(1, (j - 1 + d) / j))
}

# Convolution causale tronquée : renvoie  sum_{l=0}^{t-1} w_l x_{t-l}  pour t = 1..length(x)
# (R : convolve(x, rev(w), type = "open") est la convolution usuelle, calculée par FFT)
conv_causale <- function(x, w) {
  n <- length(x)
  convolve(x, rev(w), type = "open")[seq_len(n)]
}

# Différenciation fractionnaire (1 - B)^d x, tronquée au début de l'échantillon
diff_frac <- function(x, d) conv_causale(x, coef_diff_frac(d, length(x)))

# ---- 2. Poids MA(infini) et AR(infini) -----------------------------------------------------------------
# psi_0..psi_{H-1} du filtre Theta(B)/Phi(B)
psi_arma <- function(ar, ma, H) {
  if (H <= 1) return(1)
  if (!length(ar) && !length(ma)) return(c(1, rep(0, H - 1)))
  c(1, ARMAtoMA(ar = if (length(ar)) ar else 0, ma = if (length(ma)) ma else 0, lag.max = H - 1))
}

# Poids MA(infini) du FARIMA : X_t = sum_j psi_j eps_{t-j},  psi = [Theta/Phi] * (1-B)^(-d)
psi_farima <- function(ar, ma, d, H) conv_causale(psi_arma(ar, ma, H), coef_ma_frac(d, H))

# Poids AR(infini) : sum_j pi_j (X_{t-j} - mu) = eps_t, pi = [Phi * (1-B)^d] / Theta
# (pi_0 = 1). La prévision à un pas est  - sum_{j>=1} pi_j (X_{t+1-j} - mu).
pi_farima <- function(ar, ma, d, K) {
  poly_phi <- c(1, -ar, rep(0, K))[seq_len(K + 1)]                    # Phi(B), tronqué/complété à K+1 termes
  u <- conv_causale(coef_diff_frac(d, K + 1), poly_phi)                 # Phi(B) (1-B)^d
  if (length(ma)) as.numeric(stats::filter(u, -ma, method = "recursive")) else u   # division par Theta(B)
}

# ---- 3. Densité spectrale (forme) ----------------------------------------------------------------------------
# f(lambda) = sigma2/(2 pi) * g(lambda),   g = |Theta(e^{-i l})|^2 / |Phi(e^{-i l})|^2 * (2 sin(l/2))^(-2d)
shape_spectrale <- function(lambda, d, ar = numeric(0), ma = numeric(0)) {
  z <- exp(-1i * lambda)
  num <- 1 + (if (length(ma)) drop(outer(z, seq_along(ma), "^") %*% ma) else 0)
  den <- 1 - (if (length(ar)) drop(outer(z, seq_along(ar), "^") %*% ar) else 0)
  Mod(num)^2 / Mod(den)^2 * (2 * sin(lambda / 2))^(-2 * d)
}

# ---- 4. Autocovariances exactes -------------------------------------------------------------------------------------
# (a) Bruit intégré fractionnaire FI(d), sigma2 = 1 :
#       gamma(0) = Gamma(1-2d) / Gamma(1-d)^2 ,   gamma(h) = gamma(h-1) (h-1+d) / (h-d)
fi_acov <- function(d, H) {
  g0 <- exp(lgamma(1 - 2 * d) - 2 * lgamma(1 - d))
  i <- seq_len(H)
  g0 * cumprod(c(1, (i - 1 + d) / (i - d)))          # lags 0..H
}

# (b) Partie ARMA (sigma2 = 1) : gamma_ARMA(h) = sum_j psi_j psi_{j+h}, h = 0..M
acov_arma <- function(ar, ma, M, L = 5000) {
  psi <- psi_arma(ar, ma, L + 1)
  vapply(0:M, function(h) sum(psi[seq_len(L + 1 - h)] * psi[(1 + h):(L + 1)]), numeric(1))
}

# Horizon M au-delà duquel gamma_ARMA est négligeable (décroissance géométrique en rho^h)
horizon_acov_arma <- function(ar, ma, tol = 1e-12, M_max = 3000) {
  q <- length(ma)
  if (!length(ar)) return(q + 1L)
  rts <- polyroot(c(1, -ar))
  if (!length(rts)) return(q + 1L)
  rho <- max(Mod(1 / rts))
  if (rho < 1e-6) return(q + 1L)
  M <- ceiling(log(tol) / log(rho)) + q
  as.integer(max(q + 1L, min(M, M_max)))
}

# (c) FARIMA(p,d,q) : X = [Theta/Phi] U avec U = FI(d)  ==>  gamma_X(k) = sum_h gamma_ARMA(h) gamma_U(k - h)
#     Formule EXACTE (aucune approximation hors troncature géométrique de gamma_ARMA).
acov_farima <- function(n, d, ar = numeric(0), ma = numeric(0), sigma2 = 1) {
  if (!length(ar) && !length(ma)) return(sigma2 * fi_acov(d, n - 1))
  M  <- horizon_acov_arma(ar, ma)
  a  <- acov_arma(ar, ma, M)
  gu <- fi_acov(d, n - 1 + M)
  hh <- -M:M; ah <- a[abs(hh) + 1]
  g <- vapply(0:(n - 1), function(k) sum(ah * gu[abs(k - hh) + 1]), numeric(1))
  sigma2 * g
}

acf_farima <- function(n, d, ar = numeric(0), ma = numeric(0)) {
  g <- acov_farima(n, d, ar, ma); g / g[1]
}

# ---- 5. Durbin-Levinson : innovations et variances d'innovation ----------------------------------------------------------
# Pour un vecteur x ~ N(0, Gamma), calcule les erreurs de prévision à un pas e_t = x_t - x_hat_t
# et leurs variances v_t (relatives à sigma2). Le vecteur de uns passe dans la même récursion
# (e1) : cela donne gratuitement la moyenne GLS  mu = (1' G^-1 x)/(1' G^-1 1) = sum(e1 e / v)/sum(e1^2 / v).
# Complexité O(n^2).
dl_innovations <- function(x, gam) {
  n <- length(x)
  v <- numeric(n); e <- numeric(n); e1 <- numeric(n)
  v[1] <- gam[1]; e[1] <- x[1]; e1[1] <- 1
  phi <- numeric(0)
  for (k in seq_len(n - 1)) {
    acc <- if (k > 1) sum(phi * gam[k:2]) else 0
    kap <- (gam[k + 1] - acc) / v[k]                  # corrélation partielle d'ordre k
    phi <- c(phi - kap * rev(phi), kap)
    v[k + 1]  <- v[k] * (1 - kap^2)
    e[k + 1]  <- x[k + 1] - sum(phi * x[k:1])
    e1[k + 1] <- 1 - sum(phi)
  }
  list(e = e, e1 = e1, v = v)
}

# ---- 6. Simulation --------------------------------------------------------------------------------------------------------
# Davies-Harte : simulation gaussienne EXACTE par plongement circulant de la matrice de covariance.
sim_farima <- function(n, d, ar = numeric(0), ma = numeric(0), sigma = 1, mu = 0) {
  gam <- sigma^2 * acov_farima(n + 1, d, ar, ma)               # lags 0..n
  c_ <- c(gam[1:(n + 1)], gam[n:2])                           # première ligne du circulant (longueur 2n)
  lam <- Re(fft(c_))
  if (min(lam) < -1e-8 * max(lam)) {
    warning("Plongement circulant non défini positif : repli sur MA(infini) tronqué.")
    N <- n + 2 * n
    return(mu + sim_farima_innov(rnorm(N, sd = sigma), d, ar, ma)[(2 * n + 1):N])
  }
  m <- length(c_)
  w <- complex(real = rnorm(m), imaginary = rnorm(m))
  y <- fft(sqrt(pmax(lam, 0)) * w) / sqrt(m)
  mu + Re(y)[seq_len(n)]
}

# Filtrage d'innovations données (utile pour le bootstrap) : X_t = sum_{j<t} psi_j innov_{t-j}.
# Troncature du MA(infini) : on jette un "burn-in" au début côté appelant.
sim_farima_innov <- function(innov, d, ar = numeric(0), ma = numeric(0)) {
  conv_causale(innov, psi_farima(ar, ma, d, length(innov)))
}

# ---- 7. Reparamétrisation (optimisation sans contrainte) -------------------------------------------------------------------------
# AR : corrélations partielles r_k in (-1,1)  ->  coefficients (stationnarité garantie)
pacf_vers_ar <- function(r) {
  p <- length(r); if (p == 0) return(numeric(0))
  phi <- r[1]
  if (p > 1) for (k in 2:p) phi <- c(phi - r[k] * rev(phi), r[k])
  phi
}
ar_vers_pacf <- function(phi) {
  p <- length(phi); if (p == 0) return(numeric(0))
  r <- numeric(p); a <- phi
  for (k in p:1) {
    r[k] <- a[k]
    if (k > 1) a <- (a[1:(k - 1)] + r[k] * rev(a[1:(k - 1)])) / (1 - r[k]^2)
  }
  r
}

# u (R^k) -> (d, ar, ma).  d = d_max * tanh(u_d) ;  ar = pacf_vers_ar(tanh(u_ar)) ;  ma = -pacf_vers_ar(tanh(u_ma))
# (Theta(B) = 1 + sum theta_j B^j  est inversible  <=>  1 - sum(-theta_j) B^j  a ses racines hors du disque).
u_vers_param <- function(u, p, q, d_fixe = NULL) {
  k <- 0
  if (is.null(d_fixe)) { d <- CONFIG$d_max * tanh(u[1]); k <- 1 } else d <- d_fixe
  list(d  = d,
       ar = if (p > 0) pacf_vers_ar(tanh(u[k + seq_len(p)])) else numeric(0),
       ma = if (q > 0) -pacf_vers_ar(tanh(u[k + p + seq_len(q)])) else numeric(0))
}

param_vers_u <- function(d, ar, ma, d_fixe = NULL) {
  cl <- function(z) pmax(pmin(z, 0.999), -0.999)
  c(if (is.null(d_fixe)) atanh(cl(d / CONFIG$d_max)),
    if (length(ar)) atanh(cl(ar_vers_pacf(ar))),
    if (length(ma)) atanh(cl(ar_vers_pacf(-ma))))
}

# Vecteur "naturel" eta = (d, ar, ma) <-> composantes
eta_vers_param <- function(eta, p, q, d_fixe = NULL) {
  k <- 0
  if (is.null(d_fixe)) { d <- eta[1]; k <- 1 } else d <- d_fixe
  list(d = d, ar = if (p > 0) eta[k + seq_len(p)] else numeric(0),
       ma = if (q > 0) eta[k + p + seq_len(q)] else numeric(0))
}

# Racines : modules des racines de Phi et de Theta (>1 requis) et distance minimale AR/MA
# (une quasi-annulation de racines signale un sur-paramétrage ARMA).
verifier_racines <- function(ar, ma) {
  rac_ar <- if (length(ar)) Mod(polyroot(c(1, -ar))) else numeric(0)
  rac_ma <- if (length(ma)) Mod(polyroot(c(1, ma)))  else numeric(0)
  dist <- NA_real_
  if (length(ar) && length(ma))
    dist <- min(abs(outer(polyroot(c(1, -ar)), polyroot(c(1, ma)), "-")))
  list(module_racines_ar = rac_ar, module_racines_ma = rac_ma,
       stationnaire = all(rac_ar > 1), inversible = all(rac_ma > 1),
       distance_min_ar_ma = dist)
}
