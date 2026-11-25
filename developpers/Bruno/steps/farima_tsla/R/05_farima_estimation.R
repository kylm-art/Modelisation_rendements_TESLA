# =============================================================================
# 05_farima_estimation.R  --  Estimation GAUSSIENNE paramétrique du FARIMA(p, d, q)
# =============================================================================
# Deux estimateurs gaussiens (quasi-vraisemblance, donc consistants même si le bruit n'est pas
# gaussien, à condition que le bruit soit un bon bruit blanc) :
#
#  (A) WHITTLE  : vraisemblance gaussienne approchée dans le domaine fréquentiel.
#        Q(eta) = log( mean_j I_j / g_j(eta) ) + mean_j log g_j(eta),  j = 1..floor((n-1)/2)
#      sigma2 est concentré : sigma2_hat = mean_j I_j / g_j. Rapide (FFT), indépendant de la moyenne.
#
#  (B) MLE EXACT : vraisemblance gaussienne exacte via Durbin-Levinson, moyenne mu profilée par GLS,
#      sigma2 concentré. Coût O(n^2) par évaluation : à réserver aux modèles présélectionnés.
#
# Écarts-types asymptotiques (communs aux deux) :  Cov(eta_hat) = 2 D^{-1} / n,
#      D = (1/m) sum_j  grad log g_j  grad log g_j'      (matrice d'information de Whittle).
#   Cas particulier FI(d) : Var(d_hat) = 6 / (pi^2 n).
# Ces écarts-types supposent un bruit blanc FORT (i.i.d.) ; en présence d'effets ARCH ils sont
# trop optimistes -> voir bootstrap par blocs (module 10).
# =============================================================================

# ---- Utilitaires d'optimisation ----------------------------------------------------------------------------
optimiser_multistart <- function(fn, starts) {
  best <- NULL
  for (s in starts) {
    o <- tryCatch(optim(s, fn, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) {   # repli : Nelder-Mead (robuste, plus lent)
    for (s in starts) {
      o <- tryCatch(optim(s, fn, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10)),
                    error = function(e) NULL)
      if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
    }
  }
  best
}

construire_starts <- function(p, q, d0, d_fixe) {
  zeros <- rep(0, p + q)
  du <- function(d) if (is.null(d_fixe)) atanh(max(min(d / CONFIG$d_max, 0.99), -0.99)) else NULL
  list(c(du(d0), zeros), c(du(0), zeros), c(du(0.1), rep(0.3, p + q)))
}

# ---- Écarts-types asymptotiques (matrice D de Whittle) ------------------------------------------------------------
cov_asymptotique <- function(eta, n, p, q, d_fixe = NULL, m = floor((n - 1) / 2), h = 1e-5) {
  k <- length(eta); if (k == 0) return(NULL)
  lam <- 2 * pi * seq_len(m) / n
  logg <- function(e) { pr <- eta_vers_param(e, p, q, d_fixe)
    log(shape_spectrale(lam, pr$d, pr$ar, pr$ma)) }
  J <- vapply(seq_len(k), function(i) {
    ep <- em <- eta; ep[i] <- ep[i] + h; em[i] <- em[i] - h
    (logg(ep) - logg(em)) / (2 * h) }, numeric(m))
  D <- crossprod(J) / m
  tryCatch(2 * solve(D) / n, error = function(e) matrix(NA_real_, k, k))
}

# ---- Objet "farima_fit" ----------------------------------------------------------------------------------------------------
nouveau_fit <- function(methode, x, p, q, d_fixe, pr, mu, sigma2, V, converge, crit) {
  eta <- c(if (is.null(d_fixe)) pr$d, pr$ar, pr$ma)
  noms <- c(if (is.null(d_fixe)) "d", if (p > 0) paste0("ar", seq_len(p)), if (q > 0) paste0("ma", seq_len(q)))
  noms <- noms %||% character(0)
  se <-if (!is.null(V) && length(eta)) sqrt(pmax(diag(V), 0)) else rep(NA_real_, length(eta))
  coefs <- data.frame(param = noms, estimation = eta, se = se, z = eta / se,
                      p_value = 2 * pnorm(-abs(eta / se)), stringsAsFactors = FALSE)
  structure(list(methode = methode, p = p, q = q, d = pr$d, ar = pr$ar, ma = pr$ma, d_fixe = d_fixe,
                 mu = mu, sigma2 = sigma2, n = length(x), coefs = coefs, cov = V,
                 criteres = crit, converge = converge, x = x), class = "farima_fit")
}

print.farima_fit <- function(x, digits = 4, ...) {
  cat(sprintf("FARIMA(%d, %s, %d)  --  estimation %s  (n = %d)\n", x$p,
              if (is.null(x$d_fixe)) "d" else format(x$d_fixe), x$q, x$methode, x$n))
  print(arrondir_df(x$coefs, digits), row.names = FALSE)
  cat(sprintf("mu = %.6f   sigma2 = %.6g   convergence : %s\n", x$mu, x$sigma2,
              ifelse(isTRUE(x$converge), "oui", "NON")))
  print(round(x$criteres, 3))
  invisible(x)
}

coef.farima_fit <- function(object, ...) setNames(object$coefs$estimation, object$coefs$param)

# =============================================================================
# (A) WHITTLE
# =============================================================================
obj_whittle <- function(u, prep, p, q, d_fixe = NULL) {
  pr <- u_vers_param(u, p, q, d_fixe)
  g <- shape_spectrale(prep$lam, pr$d, pr$ar, pr$ma)
  if (any(!is.finite(g)) || any(g <= 0)) return(1e10)
  val <- log(mean(prep$I / g)) + mean(log(g))
  if (is.finite(val)) val else 1e10
}

estim_whittle <- function(x, p = 0, q = 0, d_fixe = NULL, d_init = NULL) {
  n <- length(x); per <- periodogramme(x)
  prep <- list(lam = per$lambda, I = per$I)
  d0 <- 0.1
  if (is.null(d_fixe)) {
    d0 <- d_init %||% tryCatch(est_local_whittle(x, 0.65)$d, error = function(e) 0.1)
    d0 <- max(min(d0, 0.4), -0.4)
  }
  fn <- function(u) obj_whittle(u, prep, p, q, d_fixe)

  if (is.null(d_fixe) || p + q > 0) {
    opt <- optimiser_multistart(fn, construire_starts(p, q, d0, d_fixe))
    if (is.null(opt)) stop("Whittle : l'optimisation a échoué (p = ", p, ", q = ", q, ").", call. = FALSE)
    u_hat <- opt$par; converge <- opt$convergence == 0
  } else { u_hat <- numeric(0); converge <- TRUE }          # bruit blanc pur : rien à estimer

  pr <- u_vers_param(u_hat, p, q, d_fixe)
  g  <- shape_spectrale(prep$lam, pr$d, pr$ar, pr$ma)
  s2 <- mean(prep$I / g)
  k  <- length(u_hat)
  eta <- c(if (is.null(d_fixe)) pr$d, pr$ar, pr$ma)
  V <- if (k > 0) cov_asymptotique(eta, n, p, q, d_fixe) else NULL
  crit <- c(sigma2 = s2,
            AIC = n * log(s2) + 2 * k,
            BIC = n * log(s2) + k * log(n),
            HQ  = n * log(s2) + 2 * k * log(log(n)),
            k = k)
  nouveau_fit("Whittle", x, p, q, d_fixe, pr, mu = mean(x), sigma2 = s2, V = V, converge = converge, crit = crit)
}

# =============================================================================
# (B) MAXIMUM DE VRAISEMBLANCE GAUSSIEN EXACT
# =============================================================================
# -2 log-vraisemblance concentrée (hors constante n log 2pi + n) :  n log sigma2_hat + sum log v_t
nll_exact_param <- function(pr, x) {
  n <- length(x)
  gam <- acov_farima(n, pr$d, pr$ar, pr$ma)
  if (any(!is.finite(gam)) || gam[1] <= 0) return(list(val = 1e10))
  inn <- dl_innovations(x, gam)
  if (any(!is.finite(inn$v)) || any(inn$v <= 0)) return(list(val = 1e10))
  mu <- sum(inn$e1 * inn$e / inn$v) / sum(inn$e1^2 / inn$v)       # moyenne GLS (profilée)
  e  <- inn$e - mu * inn$e1
  s2 <- mean(e^2 / inn$v)
  val <- n * log(s2) + sum(log(inn$v))
  list(val = if (is.finite(val)) val else 1e10, mu = mu, sigma2 = s2)
}

estim_mle_exact <- function(x, p = 0, q = 0, d_fixe = NULL, init = NULL, hessien = FALSE) {
  n <- length(x)
  w <- init %||% estim_whittle(x, p, q, d_fixe)               # point de départ : Whittle
  u0 <- param_vers_u(w$d, w$ar, w$ma, d_fixe)
  fn <- function(u) nll_exact_param(u_vers_param(u, p, q, d_fixe), x)$val

  if (length(u0) > 0) {
    opt <- optimiser_multistart(fn, list(u0))
    if (is.null(opt)) stop("MLE exact : l'optimisation a échoué.", call. = FALSE)
    u_hat <- opt$par; converge <- opt$convergence == 0
  } else { u_hat <- numeric(0); converge <- TRUE }

  pr  <- u_vers_param(u_hat, p, q, d_fixe)
  fin <- nll_exact_param(pr, x)
  k <- length(u_hat); eta <- c(if (is.null(d_fixe)) pr$d, pr$ar, pr$ma)
  V <- if (k > 0) cov_asymptotique(eta, n, p, q, d_fixe) else NULL
  ll <- -0.5 * (n * log(2 * pi) + fin$val + n)
  kt <- k + 2                                                  # + mu + sigma2
  crit <- c(sigma2 = fin$sigma2, loglik = ll, AIC = -2 * ll + 2 * kt, BIC = -2 * ll + kt * log(n),
            HQ = -2 * ll + 2 * kt * log(log(n)), k = k)
  fit <- nouveau_fit("MLE exact", x, p, q, d_fixe, pr, mu = fin$mu, sigma2 = fin$sigma2, V = V,
                     converge = converge, crit = crit)
  if (hessien && k > 0) fit$cov_hessien <- cov_hessienne(fit)
  fit
}

# Covariance par hessienne numérique de la log-vraisemblance exacte (contrôle indépendant de la formule
# asymptotique). Cov = (Hessienne de -loglik)^-1 = 2 * (Hessienne de n log s2 + sum log v)^-1.
cov_hessienne <- function(fit) {
  eta <- c(if (is.null(fit$d_fixe)) fit$d, fit$ar, fit$ma); k <- length(eta)
  fn <- function(e) { pr <- eta_vers_param(e, fit$p, fit$q, fit$d_fixe); nll_exact_param(pr, fit$x)$val }
  H <- tryCatch(optimHess(eta, fn, control = list(ndeps = rep(1e-4, k))), error = function(e) NULL)
  if (is.null(H)) return(matrix(NA_real_, k, k))
  tryCatch(2 * solve(H), error = function(e) matrix(NA_real_, k, k))
}

# ---- Point d'entrée unique ------------------------------------------------------------------------------------------------------
estim_farima <- function(x, p = 0, q = 0, methode = c("whittle", "exact"), d_fixe = NULL, ...) {
  methode <- match.arg(methode)
  if (methode == "whittle") estim_whittle(x, p, q, d_fixe, ...) else estim_mle_exact(x, p, q, d_fixe, ...)
}

# =============================================================================
# Tests sur d
# =============================================================================
# (1) Wald : d = 0 avec l'écart-type asymptotique. (2) Rapport de vraisemblance exact : FARIMA(p,d,q) contre ARMA(p,q).
# Sous H0 : d = 0 (point intérieur de (-0.5, 0.5)) LR ~ chi2(1).
test_d_nul <- function(x, p = 0, q = 0, fit_libre = NULL, fit_d0 = NULL) {
  f1 <- fit_libre %||% estim_mle_exact(x, p, q)
  f0 <- fit_d0 %||% estim_mle_exact(x, p, q, d_fixe = 0)
  wald_z <- f1$coefs$estimation[f1$coefs$param == "d"] / f1$coefs$se[f1$coefs$param == "d"]
  lr <- 2 * (f1$criteres["loglik"] - f0$criteres["loglik"])
  data.frame(test = c("Wald (d = 0)", "Rapport de vraisemblance (d = 0)"),
             statistique = c(wald_z^2, unname(lr)), df = 1,
             p_value = c(pchisq(wald_z^2, 1, lower.tail = FALSE), pchisq(max(lr, 0), 1, lower.tail = FALSE)),
             stringsAsFactors = FALSE)
}

# Erreur-type de la MOYENNE sous le modèle ajusté (la moyenne converge à la vitesse n^(d - 1/2) :
# l'écart-type usuel sd/sqrt(n) est FAUX sous mémoire longue).
se_moyenne_modele <- function(fit) {
  n <- fit$n; g <- acov_farima(n, fit$d, fit$ar, fit$ma, sigma2 = fit$sigma2)
  h <- 1:(n - 1)
  sqrt((g[1] + 2 * sum((1 - h / n) * g[h + 1])) / n)
}

# =============================================================================
# Innovations (résidus) EXACTES : e_t / sqrt(v_t) avec les innovations de Durbin-Levinson
# =============================================================================
residus_farima <- function(fit) {
  gam <- acov_farima(fit$n, fit$d, fit$ar, fit$ma)
  inn <- dl_innovations(fit$x - fit$mu, gam)
  eps <- inn$e / sqrt(inn$v)                      # estimation de eps_t (variance ~ sigma2)
  list(eps = eps, z = eps / sqrt(fit$sigma2), v = inn$v)
}

# =============================================================================
# Contrôles croisés
# =============================================================================
comparer_estimateurs <- function(x, p = 0, q = 0, inclure_exact = TRUE) {
  lignes <- list()
  w <- estim_whittle(x, p, q)
  lignes$whittle <- c(coef(w)["d"], coef(w)[-1])
  if (inclure_exact) { e <- estim_mle_exact(x, p, q, init = w); lignes$exact <- c(coef(e)["d"], coef(e)[-1]) }
  if (pkg_dispo("fracdiff")) {
    f <- tryCatch(fracdiff::fracdiff(x - mean(x), nar = p, nma = q), error = function(e) NULL)
    if (!is.null(f)) lignes$fracdiff <- c(d = unname(f$d), setNames(f$ar, if (p) paste0("ar", seq_len(p))),
                                        setNames(f$ma, if (q) paste0("ma", seq_len(q))))
  }
  noms <- unique(unlist(lapply(lignes, names)))
  M <- t(sapply(lignes, function(l) l[noms])); colnames(M) <- noms
  rownames(M) <- c(whittle = "Whittle", exact = "MLE exact", fracdiff = "fracdiff (Haslett-Raftery)")[names(lignes)]
  M
}
