# =============================================================================
# 10_extensions.R  --  (a) bootstrap par blocs des écarts-types, (b) FARIMA-GARCH / FIGARCH, (c) VaR
# =============================================================================
# Pourquoi des extensions ?
#  - Les écarts-types de Whittle/MLE supposent un bruit i.i.d. Pour des rendements (kurtosis élevée,
#    clusters de volatilité), le bruit est seulement non corrélé : on rééchantillonne donc les RÉSIDUS PAR
#    BLOCS (qui conservent la dépendance en variance) pour obtenir des écarts-types honnêtes.
#  - Le FARIMA modélise la MOYENNE conditionnelle ; il ne dit rien de la variance conditionnelle.
#    Le modèle complet est FARIMA-GARCH (ou FARIMA-FIGARCH si la volatilité elle-même est à mémoire longue).
# =============================================================================

# ---- (a) Bootstrap par blocs mobiles des résidus ---------------------------------------------------------------------
# Schéma : 1) résidus estimés eps_hat ; 2) rééchantillonnage par blocs de longueur L ~ n^(1/3) ;
#          3) reconstruction x* = mu + sum_j psi_j eps*_{t-j} (avec burn-in) ; 4) ré-estimation.
# Approximation : le MA(infini) est tronqué (burn-in = n) -> léger biais sur la partie très basse fréquence.
bootstrap_blocs_farima <- function(fit, B = 200, longueur_bloc = NULL, burn = NULL, graine = CONFIG$graine,
                                   progression = TRUE) {
  set.seed(graine)
  n <- fit$n; eps <- residus_farima(fit)$eps; eps <- eps - mean(eps)
  L <- longueur_bloc %||% ceiling(n^(1 / 3)); burn <- burn %||% n; N <- n + burn
  psi <- psi_farima(fit$ar, fit$ma, fit$d, N)
  est <- matrix(NA_real_, B, nrow(fit$coefs), dimnames = list(NULL, fit$coefs$param))
  for (b in seq_len(B)) {
    if (progression && b %% 50 == 0) message(sprintf("  bootstrap : %d / %d", b, B))
    nb <- ceiling(N / L); deb <- sample.int(n - L + 1, nb, replace = TRUE)
    e_star <- unlist(lapply(deb, function(s) eps[s:(s + L - 1)]))[seq_len(N)]
    xs <- fit$mu + conv_causale(e_star, psi)[(burn + 1):N]
    f <- tryCatch(estim_whittle(xs, fit$p, fit$q, fit$d_fixe), error = function(e) NULL)
    if (!is.null(f)) est[b, ] <- f$coefs$estimation
  }
  tab <- data.frame(param = fit$coefs$param, estimation = fit$coefs$estimation,
                    se_asymptotique = fit$coefs$se,
                    se_bootstrap = apply(est, 2, sd, na.rm = TRUE),
                    ic95_bas = apply(est, 2, quantile, 0.025, na.rm = TRUE),
                    ic95_haut = apply(est, 2, quantile, 0.975, na.rm = TRUE),
                    reussites = colSums(!is.na(est)), row.names = NULL)
  tab$ratio_se <- tab$se_bootstrap / tab$se_asymptotique
  list(estim = est, tab = tab, longueur_bloc = L)
}

plot_bootstrap_d <- function(bs, d_hat = NULL) {
  par_std()
  e <- bs$estim[, "d"]; e <- e[is.finite(e)]
  hist(e, breaks = 30, col = "grey85", border = "white", freq = FALSE, xlab = "d*", main = "Distribution bootstrap de d")
  abline(v = d_hat %||% mean(e), col = COULEURS["rouge"], lwd = 2)
  abline(v = 0, lty = 2)
}

# ---- (b) FARIMA-GARCH / FARIMA-FIGARCH via rugarch -----------------------------------------------------------------------------
# rugarch : arfima = TRUE ajoute la différenciation fractionnaire (1-B)^d à la moyenne ARMA.
# variance = "sGARCH" (mémoire courte), "gjrGARCH" (asymétrie : effet de levier),
#            "fiGARCH" (mémoire longue dans la variance).
ajuster_farima_garch <- function(y, p = 0, q = 0, ordre_garch = c(1, 1),
                                 variance = c("sGARCH", "gjrGARCH", "fiGARCH"),
                                 distribution = "std", arfima = TRUE) {
  if (!pkg_dispo("rugarch")) { message("Package 'rugarch' non installé : étape ignorée."); return(NULL) }
  variance <- match.arg(variance)
  spec <- rugarch::ugarchspec(
    variance.model = list(model = variance, garchOrder = ordre_garch),
    mean.model = list(armaOrder = c(p, q), include.mean = TRUE, arfima = arfima),
    distribution.model = distribution)
  fit <- tryCatch(rugarch::ugarchfit(spec, data = y, solver = "hybrid"),
                  error = function(e) { message("ugarchfit a échoué : ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  z <- as.numeric(rugarch::residuals(fit, standardize = TRUE))
  list(fit = fit, spec = spec, z = z, variance = variance, distribution = distribution,
       coefs = fit@fit$robust.matcoef,                     # écarts-types robustes (sandwich)
       ic = rugarch::infocriteria(fit)[, 1],
       sigma = as.numeric(rugarch::sigma(fit)))
}

diagnostics_garch <- function(res, lags = c(10, 20)) {
  z <- res$z
  list(portmanteau_z   = test_portmanteau(z, lags, fitdf = 0),
       portmanteau_z2  = do.call(rbind, lapply(lags, function(m) {
         t <- Box.test(z^2, lag = m, type = "Ljung-Box"); data.frame(retards = m, stat = unname(t$statistic), p = t$p.value) })),
       arch = arch_lm(z, 10),
       normalite = data.frame(asymetrie = asymetrie(z), kurtosis = kurtosis(z)))
}

# Comparaison de plusieurs variantes de variance (même moyenne FARIMA)
comparer_variances <- function(y, p = 0, q = 0, variances = c("sGARCH", "gjrGARCH", "fiGARCH"),
                               distribution = "std") {
  fits <- lapply(variances, function(v) ajuster_farima_garch(y, p, q, c(1, 1), v, distribution))
  names(fits) <- variances
  ok <- !vapply(fits, is.null, logical(1))
  tab <- do.call(rbind, lapply(names(fits)[ok], function(v) data.frame(variance = v, t(fits[[v]]$ic))))
  list(fits = fits, table = tab)
}

# ---- (c) Prévision conditionnelle et Value-at-Risk ------------------------------------------------------------------------------------------
var_farima_garch <- function(res, h = 10, niveaux = c(0.01, 0.05)) {
  fc  <- rugarch::ugarchforecast(res$fit, n.ahead = h)
  mu  <- as.numeric(rugarch::fitted(fc)); sg <- as.numeric(rugarch::sigma(fc))
  cf  <- rugarch::coef(res$fit); shape <- if ("shape" %in% names(cf)) cf[["shape"]] else NULL
  out <- data.frame(h = seq_len(h), moyenne = mu, volatilite_jour = sg, volatilite_annuelle = sg * sqrt(CONFIG$jours_an))
  for (a in niveaux) {
    qa <- rugarch::qdist(res$distribution, p = a, mu = 0, sigma = 1, shape = shape %||% 5)
    out[[paste0("VaR_", (1 - a) * 100)]] <- mu + sg * qa
  }
  out
}

# Test de Kupiec (couverture inconditionnelle) sur des violations de VaR
kupiec <- function(violations, alpha) {
  T <- length(violations); x <- sum(violations); pi_hat <- x / T
  ll0 <- (T - x) * log(1 - alpha) + x * log(alpha)
  ll1 <- ifelse(x == 0 || x == T, 0, (T - x) * log(1 - pi_hat) + x * log(pi_hat))
  lr <- -2 * (ll0 - ll1)
  data.frame(taux_observe = pi_hat, taux_theorique = alpha, violations = x, LR = lr,
             p_value = pchisq(lr, 1, lower.tail = FALSE))
}

backtest_var_insample <- function(res, y, niveaux = c(0.01, 0.05)) {
  do.call(rbind, lapply(niveaux, function(a) {
    cf <- rugarch::coef(res$fit); shape <- if ("shape" %in% names(cf)) cf[["shape"]] else NULL
    qa <- rugarch::qdist(res$distribution, p = a, mu = 0, sigma = 1, shape = shape %||% 5)
    mu_t <- as.numeric(rugarch::fitted(res$fit))
    var_t <- mu_t + res$sigma * qa
    cbind(niveau_VaR = 1 - a, kupiec(y < var_t, a))
  }))
}
