# =============================================================================
# 08_prevision.R  --  Prévision à horizon h avec intervalles, et vérification par le prédicteur exact
# =============================================================================
# Prévision à partir de la représentation AR(infini) :
#     sum_{j>=0} pi_j (X_{t-j} - mu) = eps_t   (pi_0 = 1)
#   =>  X_hat_{t+h} = mu - sum_{j>=1} pi_j (X_hat_{t+h-j} - mu)   (récursif, avec X_hat = X pour le passé observé)
# Erreur quadratique moyenne à l'horizon h (poids MA(infini) psi_j) :
#     MSE_h = sigma2 * sum_{j=0}^{h-1} psi_j^2
# Sous mémoire longue, pi_j ~ j^(-d-1) décroît LENTEMENT : il faut remonter loin dans le passé (K grand).
# Les intervalles ci-dessous sont gaussiens ET INCONDITIONNELS : ils ne s'adaptent pas à la volatilité
# du moment (effet ARCH) ; pour des intervalles conditionnels, voir FARIMA-GARCH (module 10).
# =============================================================================

prevoir_farima <- function(fit, h = 20, niveaux = c(0.80, 0.95), K = NULL) {
  n <- fit$n; K <- min(K %||% n, n)
  xc  <- fit$x - fit$mu
  piw <- pi_farima(fit$ar, fit$ma, fit$d, K)                 # pi_0..pi_K
  y <- c(tail(xc, K), numeric(h))
  for (i in seq_len(h)) {
    idx <- K + i
    y[idx] <- -sum(piw[2:(K + 1)] * y[(idx - 1):(idx - K)])
  }
  pred <- fit$mu + y[K + seq_len(h)]
  psi  <- psi_farima(fit$ar, fit$ma, fit$d, h)
  se   <- sqrt(fit$sigma2 * cumsum(psi^2))
  out <- data.frame(h = seq_len(h), prevision = pred, se = se)
  for (nv in niveaux) {
    z <- qnorm(0.5 + nv / 2)
    out[[paste0("bas_", nv * 100)]]  <- pred - z * se
    out[[paste0("haut_", nv * 100)]] <- pred + z * se
  }
  out
}

# Prévision du PRIX : P_{T+h} = P_T exp( r_{T+1} + ... + r_{T+h} )
# Variance de la somme cumulée : MSE_cum(h) = sigma2 * sum_{m=0}^{h-1} Psi_m^2,  Psi_m = sum_{l<=m} psi_l.
prevoir_prix <- function(fit, P0, h = 60, niveaux = c(0.80, 0.95), K = NULL) {
  pf  <- prevoir_farima(fit, h, niveaux, K)
  psi <- psi_farima(fit$ar, fit$ma, fit$d, h)
  se_cum <- sqrt(fit$sigma2 * cumsum(cumsum(psi)^2))
  cum <- cumsum(pf$prevision)
  out <- data.frame(h = pf$h, prix_central = P0 * exp(cum), se_log = se_cum)
  for (nv in niveaux) {
    z <- qnorm(0.5 + nv / 2)
    out[[paste0("bas_", nv * 100)]]  <- P0 * exp(cum - z * se_cum)
    out[[paste0("haut_", nv * 100)]] <- P0 * exp(cum + z * se_cum)
  }
  out
}

# ---- Prédicteur linéaire optimal EXACT (contrôle) -----------------------------------------------------------------------
# x_hat_{n+s} = c_s' Gamma_K^{-1} x_{derniers K}, avec Gamma_K de Toeplitz construite à partir des autocovariances EXACTES.
prevoir_exact <- function(fit, h = 10, K = 500) {
  n <- fit$n; K <- min(K, n)
  gam <- acov_farima(K + h, fit$d, fit$ar, fit$ma, sigma2 = fit$sigma2)
  G <- toeplitz(gam[seq_len(K)])
  xK <- tail(fit$x - fit$mu, K)
  w <- solve(G, xK)
  fit$mu + vapply(seq_len(h), function(s) sum(gam[(K + s):(s + 1)] * w), numeric(1))
}

# Compare prédicteur AR(infini) tronqué et prédicteur exact (écart max, en unités de l'écart-type des données)
verifier_prevision <- function(fit, h = 10, K = 500) {
  a <- prevoir_farima(fit, h, K = K)$prevision; e <- prevoir_exact(fit, h, K)
  data.frame(h = seq_len(h), ar_infini = a, exact = e, ecart = a - e,
             ecart_rel_sd = (a - e) / sd(fit$x))
}

# ---- Graphiques -----------------------------------------------------------------------------------------------------------------------
plot_prevision <- function(fit, pf, n_hist = 120, dates = NULL, niveau = 95) {
  n <- fit$n; idx_h <- (n - n_hist + 1):n; h <- nrow(pf)
  xs <- c(idx_h, n + seq_len(h)); lo <- pf[[paste0("bas_", niveau)]]; hi <- pf[[paste0("haut_", niveau)]]
  par_std()
  plot(xs, c(fit$x[idx_h], pf$prevision), type = "n", ylim = range(c(fit$x[idx_h], lo, hi)),
       xlab = "indice temporel", ylab = "rendement", main = paste0("Prévision FARIMA, h = ", h, " (IC ", niveau, " %)"))
  lines(idx_h, fit$x[idx_h], col = COULEURS["gris"])
  polygon(c(n + seq_len(h), rev(n + seq_len(h))), c(lo, rev(hi)), col = adjustcolor(COULEURS["bleu"], 0.2), border = NA)
  lines(n + seq_len(h), pf$prevision, col = COULEURS["bleu"], lwd = 2)
  abline(h = fit$mu, lty = 2, col = COULEURS["rouge"])
  legend("topleft", c("observé", "prévision", "moyenne estimée"), col = c(COULEURS["gris"], COULEURS["bleu"], COULEURS["rouge"]),
         lty = c(1, 1, 2), lwd = c(1, 2, 1), bty = "n", cex = 0.8)
}

plot_fan_prix <- function(prix_hist, pp, n_hist = 250) {
  p <- as.numeric(prix_hist); n <- length(p); h <- nrow(pp); idx_h <- (n - n_hist + 1):n
  par_std()
  plot(c(idx_h, n + seq_len(h)), c(p[idx_h], pp$prix_central), type = "n",
       ylim = range(c(p[idx_h], pp$bas_95, pp$haut_95)), xlab = "jours", ylab = "prix",
       main = paste0("Prix : trajectoire centrale et intervalles (h = ", h, " jours)"))
  lines(idx_h, p[idx_h], col = COULEURS["gris"], lwd = 1.5)
  polygon(c(n + pp$h, rev(n + pp$h)), c(pp$bas_95, rev(pp$haut_95)), col = adjustcolor(COULEURS["bleu"], 0.15), border = NA)
  polygon(c(n + pp$h, rev(n + pp$h)), c(pp$bas_80, rev(pp$haut_80)), col = adjustcolor(COULEURS["bleu"], 0.25), border = NA)
  lines(n + pp$h, pp$prix_central, col = COULEURS["bleu"], lwd = 2)
  legend("topleft", c("historique", "médiane (exp. du cumul)", "IC 80 %", "IC 95 %"),
         col = c(COULEURS["gris"], COULEURS["bleu"], adjustcolor(COULEURS["bleu"], 0.25), adjustcolor(COULEURS["bleu"], 0.15)),
         lty = c(1, 1, NA, NA), pch = c(NA, NA, 15, 15), bty = "n", cex = 0.8)
}
