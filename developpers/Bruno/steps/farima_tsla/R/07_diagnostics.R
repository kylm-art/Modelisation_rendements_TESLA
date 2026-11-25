# =============================================================================
# 07_diagnostics.R  --  Validation du modèle FARIMA ajusté
# =============================================================================
# Un modèle n'est utilisable pour la prévision que si ses résidus ressemblent à un bruit blanc.
# Deux notions à distinguer (et c'est crucial pour les rendements financiers) :
#   - bruit blanc FORT  : i.i.d.                       -> tests classiques (Ljung-Box) valides
#   - bruit blanc FAIBLE: seulement non corrélé (ex. GARCH : eps_t non corrélé mais eps_t^2 corrélé)
#                          -> Ljung-Box classique REJETTE TROP SOUVENT ; il faut les versions robustes
#                             (Romano-Thombs pour l'ACF, Lobato pour le portmanteau).
# =============================================================================

# ---- ACF avec bandes de confiance "forte" et "faible" -------------------------------------------------
acf_robuste <- function(z, lag.max = 30) {
  n <- length(z); zc <- z - mean(z); g0 <- mean(zc^2); h <- seq_len(lag.max)
  gam <- vapply(h, function(k) sum(zc[(k + 1):n] * zc[1:(n - k)]) / n, numeric(1))
  tau <- vapply(h, function(k) mean((zc[(k + 1):n] * zc[1:(n - k)])^2), numeric(1))
  data.frame(h = h, rho = gam / g0, se_fort = 1 / sqrt(n),
             se_faible = sqrt(tau / n) / g0, gam = gam, tau = tau)
}

plot_acf_robuste <- function(tab, titre = "ACF des résidus") {
  lim <- max(abs(tab$rho), 2 * tab$se_faible, 2 * tab$se_fort) * 1.1
  par_std()
  plot(tab$h, tab$rho, type = "h", lwd = 2, ylim = c(-lim, lim), xlab = "retard", ylab = "autocorrélation",
       main = titre, col = "grey30")
  abline(h = 0)
  abline(h = c(-2, 2) * tab$se_fort[1], col = COULEURS["bleu"], lty = 2, lwd = 1.5)
  lines(tab$h, 2 * tab$se_faible, col = COULEURS["rouge"], lwd = 2)
  lines(tab$h, -2 * tab$se_faible, col = COULEURS["rouge"], lwd = 2)
  legend("topright", c("bande bruit fort (±2/sqrt(n))", "bande bruit faible (±2 se robuste)"),
         col = COULEURS[c("bleu", "rouge")], lty = c(2, 1), lwd = c(1.5, 2), bty = "n", cex = 0.8)
}

# ---- Portmanteau : Ljung-Box classique et version robuste de Lobato --------------------------------------------
# Lobato et al. (2001) : Q* = n * sum_{h=1}^m gam_h^2 / tau_h  ~  chi2(m - fitdf) sous H0 (bruit faible),
# tau_h = mean( (z_t z_{t-h})^2 ).
test_portmanteau <- function(z, lags = c(10, 20, 40), fitdf = 0) {
  n <- length(z); t <- acf_robuste(z, max(lags))
  do.call(rbind, lapply(lags, function(m) {
    lb <- Box.test(z, lag = m, type = "Ljung-Box", fitdf = fitdf)
    Qr <- n * sum(t$gam[1:m]^2 / t$tau[1:m])
    data.frame(retards = m, LjungBox = unname(lb$statistic), p_LB = lb$p.value,
               Lobato = Qr, p_Lobato = pchisq(Qr, m - fitdf, lower.tail = FALSE))
  }))
}

# ---- Adéquation spectrale ------------------------------------------------------------------------------------------------------
# Si le modèle est bon, R_j = I_j / (sigma2 g_j) est (asymptotiquement) i.i.d. exponentiel(1) : on teste la moyenne,
# la loi (KS) et, surtout, le comportement aux BASSES fréquences, là où se joue la mémoire longue.
diag_spectral <- function(fit) {
  per <- periodogramme(fit$x); g <- shape_spectrale(per$lambda, fit$d, fit$ar, fit$ma)
  R <- per$I / (fit$sigma2 * g); m <- length(R); nb <- max(10, ceiling(0.05 * m))
  ks <- suppressWarnings(ks.test(R, "pexp"))
  list(R = R, lambda = per$lambda, g = g,
       tab = data.frame(
         indicateur = c("moyenne de I/(s2 g)  (attendu 1)", "moyenne sur les 5 % plus basses fréquences (attendu 1)",
                        "KS vs Exp(1) (p-value)"),
         valeur = c(mean(R), mean(R[seq_len(nb)]), ks$p.value),
         ecart_type_theorique = c(1 / sqrt(m), 1 / sqrt(nb), NA_real_)))
}

plot_diag_spectral <- function(fit) {
  ds <- diag_spectral(fit); per <- periodogramme(fit$x)
  par_std(1, 2)
  plot(log10(per$lambda), log10(per$I), pch = 20, cex = 0.3, col = "grey70", xlab = "log10 lambda",
       ylab = "log10 I", main = "Périodogramme et densité ajustée")
  lines(log10(per$lambda), log10(fit$sigma2 * ds$g), col = COULEURS["rouge"], lwd = 2)
  # moyennes par classes de fréquence (log-régulières) : la densité ajustée doit les traverser
  cl <- cut(log10(per$lambda), breaks = 20)
  mm <- tapply(per$I, cl, mean); ll <- tapply(log10(per$lambda), cl, mean)
  points(ll, log10(mm), pch = 19, col = COULEURS["bleu"])
  legend("topright", c("modèle", "moyennes par classes"), col = COULEURS[c("rouge", "bleu")],
         lty = c(1, NA), pch = c(NA, 19), bty = "n", cex = 0.8)
  plot(ds$lambda, stats::filter(ds$R, rep(1 / 25, 25), sides = 2), type = "l", col = COULEURS["bleu"], lwd = 2,
       xlab = "lambda", ylab = "I / (s2 g)  lissé", main = "Rapport périodogramme / modèle (attendu = 1)")
  abline(h = 1, col = COULEURS["rouge"], lty = 2)
  par_std()
}

# ---- ACF empirique contre ACF impliquée par le modèle -------------------------------------------------------------------------------
plot_acf_modele <- function(fit, lag.max = 100, titre = "") {
  emp <- acf(fit$x, lag.max = lag.max, plot = FALSE)$acf[-1]
  the <- acf_farima(lag.max + 1, fit$d, fit$ar, fit$ma)[-1]
  par_std()
  plot(seq_len(lag.max), emp, type = "h", lwd = 2, col = "grey55", xlab = "retard", ylab = "autocorrélation",
       main = paste("ACF empirique vs modèle", titre))
  lines(seq_len(lag.max), the, col = COULEURS["rouge"], lwd = 2)
  abline(h = 0); abline(h = c(-2, 2) / sqrt(fit$n), lty = 2, col = COULEURS["bleu"])
  legend("topright", c("empirique", "FARIMA ajusté"), col = c("grey55", COULEURS["rouge"]), lwd = 2, bty = "n", cex = 0.8)
}

# ---- Diagnostic complet ----------------------------------------------------------------------------------------------------------------------
diagnostics_residus <- function(fit, lags = c(10, 20, 40)) {
  r <- residus_farima(fit); z <- r$z
  fitdf <- fit$criteres["k"] %||% (fit$p + fit$q)
  fitdf <- as.integer(fitdf)
  port <- test_portmanteau(z, lags, fitdf = fitdf)
  mcl  <- do.call(rbind, lapply(lags, function(m) {
    t <- Box.test(z^2, lag = m, type = "Ljung-Box"); data.frame(retards = m, McLeodLi = unname(t$statistic), p = t$p.value) }))
  arch <- arch_lm(z, 10)
  jb   <- tseries::jarque.bera.test(z)
  tmoy <- t.test(r$eps)
  list(residus = r, portmanteau = port, mcleod_li = mcl, arch = arch,
       normalite = data.frame(asymetrie = asymetrie(z), kurtosis = kurtosis(z),
                              JB = unname(jb$statistic), p_JB = jb$p.value),
       moyenne = data.frame(moyenne_eps = mean(r$eps), p_value = tmoy$p.value),
       spectral = diag_spectral(fit)$tab)
}

plot_diagnostics <- function(fit, diag = NULL, lag.max = 40) {
  d <- diag %||% diagnostics_residus(fit); z <- d$residus$z; n <- length(z)
  par_std(2, 3)
  plot(z, type = "l", col = COULEURS["gris"], main = "Résidus standardisés", xlab = "t", ylab = "")
  abline(h = c(-3, 3), col = COULEURS["rouge"], lty = 2)
  t <- acf_robuste(z, lag.max)
  plot(t$h, t$rho, type = "h", lwd = 2, col = "grey30", ylim = range(c(t$rho, 2 * t$se_faible, -2 * t$se_faible)),
       main = "ACF résidus (bandes forte/faible)", xlab = "retard", ylab = ""); abline(h = 0)
  abline(h = c(-2, 2) * t$se_fort[1], col = COULEURS["bleu"], lty = 2)
  lines(t$h, 2 * t$se_faible, col = COULEURS["rouge"], lwd = 2); lines(t$h, -2 * t$se_faible, col = COULEURS["rouge"], lwd = 2)
  acf(z^2, lag.max = lag.max, main = "ACF des résidus au carré")
  qqnorm(z, main = "QQ-plot des résidus", pch = 20, cex = 0.4, col = COULEURS["gris"]); qqline(z, col = COULEURS["rouge"], lwd = 2)
  stats::cpgram(z, main = "Périodogramme cumulé des résidus")
  hist(z, breaks = 60, freq = FALSE, col = "grey85", border = "white", main = "Histogramme des résidus", xlab = "")
  curve(dnorm(x), add = TRUE, col = COULEURS["rouge"], lwd = 2)
  par_std()
}

# ---- Stabilité des paramètres (sous-échantillons glissants) -------------------------------------------------------------------------------------
parametres_glissants <- function(x, p = 0, q = 0, largeur = 1500, pas = 250, dates = NULL) {
  n <- length(x); debuts <- seq(1, n - largeur + 1, by = pas)
  res <- do.call(rbind, lapply(debuts, function(i) {
    f <- tryCatch(estim_whittle(x[i:(i + largeur - 1)], p, q), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    cf <- f$coefs
    data.frame(debut = i, fin = i + largeur - 1, param = cf$param, estimation = cf$estimation, se = cf$se)
  }))
  if (!is.null(dates)) { res$date_debut <- dates[res$debut]; res$date_fin <- dates[res$fin] }
  res
}

plot_parametres_glissants <- function(res) {
  noms <- unique(res$param); par_std(1, length(noms))
  for (nm in noms) {
    r <- res[res$param == nm, ]; mid <- (r$debut + r$fin) / 2
    plot(mid, r$estimation, type = "b", pch = 19, col = COULEURS["bleu"],
         ylim = range(c(r$estimation - 2 * r$se, r$estimation + 2 * r$se, 0), na.rm = TRUE),
         xlab = "milieu de fenêtre", ylab = nm, main = paste("Stabilité :", nm))
    segments(mid, r$estimation - 2 * r$se, mid, r$estimation + 2 * r$se, col = COULEURS["bleu"]); abline(h = 0, lty = 2)
  }
  par_std()
}

# ---- Feuille de contrôle (OK / ATTENTION / ECHEC) ---------------------------------------------------------------------------------------------------------------
verifications_modele <- function(fit, diag = NULL, alpha = 0.05) {
  d <- diag %||% diagnostics_residus(fit)
  rac <- verifier_racines(fit$ar, fit$ma)
  p20 <- d$portmanteau$p_Lobato[d$portmanteau$retards == 20]
  p20c <- d$portmanteau$p_LB[d$portmanteau$retards == 20]
  mcl20 <- d$mcleod_li$p[d$mcleod_li$retards == 20]
  se_ok <- all(is.finite(fit$coefs$se))
  stat <- function(ok, warn = FALSE) ifelse(ok, "OK", ifelse(warn, "ATTENTION", "ECHEC"))
  data.frame(
    verification = c("Optimisation convergée", "Stationnarité : d < 0.5 et racines AR hors du disque",
                     "Inversibilité : racines MA hors du disque", "Pas de facteur AR/MA quasi commun",
                     "Écarts-types finis (information inversible)",
                     "Résidus non corrélés (Lobato, 20 retards)", "Résidus non corrélés (Ljung-Box classique, 20 retards)",
                     "Pas d'effet ARCH résiduel (McLeod-Li, 20 retards)", "Résidus de moyenne nulle",
                     "Adéquation spectrale (KS sur I/(s2 g))", "Adéquation spectrale aux basses fréquences"),
    p_value_ou_valeur = c(NA, fit$d, NA, rac$distance_min_ar_ma, NA, p20, p20c, mcl20, d$moyenne$p_value,
                          d$spectral$valeur[3], d$spectral$valeur[2]),
    statut = c(stat(isTRUE(fit$converge)),
               stat(abs(fit$d) < 0.5 && rac$stationnaire),
               stat(rac$inversible),
               stat(is.na(rac$distance_min_ar_ma) || rac$distance_min_ar_ma >= 0.10, warn = TRUE),
               stat(se_ok, warn = TRUE),
               stat(p20 >= alpha),
               stat(p20c >= alpha, warn = TRUE),
               stat(mcl20 >= alpha, warn = TRUE),
               stat(d$moyenne$p_value >= alpha, warn = TRUE),
               stat(d$spectral$valeur[3] >= alpha, warn = TRUE),
               stat(abs(d$spectral$valeur[2] - 1) <= 3 * d$spectral$ecart_type_theorique[2], warn = TRUE)),
    stringsAsFactors = FALSE)
}
