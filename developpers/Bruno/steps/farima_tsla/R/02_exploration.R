# =============================================================================
# 02_exploration.R  --  Statistiques descriptives, stationnarité, effets ARCH,
#                       corrélogrammes (rendements, |r|, r^2)
# =============================================================================

# ---- Moments empiriques --------------------------------------------------------------
asymetrie <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
kurtosis  <- function(x) { x <- x - mean(x); mean(x^4) / mean(x^2)^2 }   # normale = 3

stats_descriptives <- function(y, jours_an = CONFIG$jours_an) {
  jb <- tseries::jarque.bera.test(y)
  data.frame(
    statistique = c("n", "moyenne quotidienne", "écart-type quotidien",
                    "rendement annualisé (moyenne x 252)", "volatilité annualisée",
                    "minimum", "maximum", "asymétrie (skewness)", "kurtosis (normale = 3)",
                    "Jarque-Bera (stat.)", "Jarque-Bera (p-value)"),
    valeur = c(length(y), mean(y), sd(y), mean(y) * jours_an, sd(y) * sqrt(jours_an),
               min(y), max(y), asymetrie(y), kurtosis(y),
               unname(jb$statistic), jb$p.value),
    stringsAsFactors = FALSE)
}

# ---- Tests de stationnarité (croisés) ------------------------------------------------------
# ADF  : H0 = racine unitaire    | PP : H0 = racine unitaire   | KPSS : H0 = stationnarité
# Lecture croisée : ADF/PP rejettent ET KPSS ne rejette pas  =>  stationnaire (d < 1).
# ATTENTION : ces tests ne distinguent PAS d = 0 de 0 < d < 0.5 (tous deux stationnaires) ;
# c'est le rôle des estimateurs de d (module 03).
tests_stationnarite <- function(x) {
  adf  <- suppressWarnings(tseries::adf.test(x))
  pp   <- suppressWarnings(tseries::pp.test(x))
  kp   <- suppressWarnings(tseries::kpss.test(x, null = "Level"))
  data.frame(
    test = c("ADF", "Phillips-Perron", "KPSS (niveau)"),
    H0 = c("racine unitaire", "racine unitaire", "stationnarité"),
    statistique = c(unname(adf$statistic), unname(pp$statistic), unname(kp$statistic)),
    p_value = c(adf$p.value, pp$p.value, kp$p.value),
    conclusion = c(ifelse(adf$p.value < 0.05, "stationnaire", "racine unitaire non rejetée"),
                   ifelse(pp$p.value  < 0.05, "stationnaire", "racine unitaire non rejetée"),
                   ifelse(kp$p.value  < 0.05, "NON stationnaire", "stationnarité non rejetée")),
    stringsAsFactors = FALSE)
}

# ---- Test d'effet ARCH (multiplicateur de Lagrange d'Engle) --------------------------------
# Régression de (x - mean)^2 sur ses `lags` retards ; stat = n * R^2 ~ chi2(lags) sous H0.
arch_lm <- function(x, lags = 10) {
  e2 <- (x - mean(x))^2; n <- length(e2)
  M <- embed(e2, lags + 1)
  fit <- lm(M[, 1] ~ M[, -1])
  stat <- (n - lags) * summary(fit)$r.squared
  data.frame(test = paste0("ARCH-LM (", lags, " retards)"), statistique = stat,
             df = lags, p_value = pchisq(stat, lags, lower.tail = FALSE))
}

# ---- Ljung-Box sur plusieurs retards (tableau) ----------------------------------------------
ljung_box_multi <- function(x, lags = c(5, 10, 20, 40), fitdf = 0, etiquette = "x") {
  do.call(rbind, lapply(lags, function(m) {
    t <- Box.test(x, lag = m, type = "Ljung-Box", fitdf = fitdf)
    data.frame(serie = etiquette, retards = m, statistique = unname(t$statistic),
               p_value = t$p.value)
  }))
}

# ---- Décroissance de l'autocorrélation : exponentielle ou hyperbolique ? ------------------------
# Mémoire courte : rho(h) ~ a * r^h   => log rho linéaire en h
# Mémoire longue : rho(h) ~ c * h^(2d-1) => log rho linéaire en log h  (pente = 2d - 1)
# On régresse sur les retards où l'autocorrélation est positive (log défini).
decroissance_acf <- function(x, lag.max = 250) {
  a <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1]
  h <- seq_along(a); ok <- a > 0
  if (sum(ok) < 15) return(list(ok = FALSE, message = "trop peu d'autocorrélations positives"))
  fit_hyp <- lm(log(a[ok]) ~ log(h[ok]))
  fit_exp <- lm(log(a[ok]) ~ h[ok])
  list(ok = TRUE, pente_loglog = unname(coef(fit_hyp)[2]),
       d_implicite = (unname(coef(fit_hyp)[2]) + 1) / 2,
       r2_hyperbolique = summary(fit_hyp)$r.squared,
       r2_exponentielle = summary(fit_exp)$r.squared,
       n_points = sum(ok), h = h[ok], rho = a[ok])
}

# ---- Graphiques -------------------------------------------------------------------------------------
plot_apercu <- function(jeu) {
  par_std(2, 1)
  plot(zoo::as.zoo(jeu$prix), log = "y", main = paste("Prix ajusté (échelle log) :", CONFIG$ticker),
       xlab = "", ylab = "", col = COULEURS["bleu"])
  plot(zoo::as.zoo(jeu$r), main = "Rendements logarithmiques", xlab = "", ylab = "",
       col = COULEURS["gris"])
  abline(h = 0, col = COULEURS["rouge"], lty = 2)
  par_std()
}

plot_distribution <- function(y) {
  par_std(1, 3)
  hist(y, breaks = 80, freq = FALSE, col = "grey85", border = "white",
       main = "Densité des rendements", xlab = "")
  curve(dnorm(x, mean(y), sd(y)), add = TRUE, col = COULEURS["rouge"], lwd = 2)
  lines(density(y), col = COULEURS["bleu"], lwd = 2)
  legend("topright", c("empirique", "normale"), col = COULEURS[c("bleu", "rouge")], lwd = 2, bty = "n", cex = 0.8)
  qqnorm(y, main = "QQ-plot normal", pch = 20, cex = 0.4, col = COULEURS["gris"]); qqline(y, col = COULEURS["rouge"], lwd = 2)
  # queues : probabilité d'un |z| > k, empirique vs normale
  k <- 1:6; z <- abs((y - mean(y)) / sd(y))
  emp <- sapply(k, function(kk) mean(z > kk)); the <- 2 * pnorm(-k)
  plot(k, log10(pmax(emp, 1e-5)), type = "b", pch = 19, col = COULEURS["bleu"],
       ylim = range(log10(pmax(c(emp, the), 1e-5))), xlab = "k (en écarts-types)",
       ylab = "log10 P(|Z| > k)", main = "Queues épaisses")
  lines(k, log10(the), type = "b", pch = 17, col = COULEURS["rouge"])
  legend("topright", c("empirique", "normale"), col = COULEURS[c("bleu", "rouge")], pch = c(19, 17), bty = "n", cex = 0.8)
  par_std()
}

plot_correlogrammes <- function(y, lag.max = 60) {
  par_std(2, 3)
  acf(y,   lag.max = lag.max, main = "ACF  r");        pacf(y,   lag.max = lag.max, main = "PACF  r")
  acf(abs(y), lag.max = lag.max, main = "ACF  |r|")
  acf(y^2, lag.max = lag.max, main = "ACF  r^2")
  pacf(abs(y), lag.max = lag.max, main = "PACF  |r|")
  # ACF longue portée de |r| : lente décroissance = signature de mémoire longue
  acf(abs(y), lag.max = 400, main = "ACF  |r|  (400 retards)")
  par_std()
}

plot_decroissance <- function(dec, titre = "") {
  if (!isTRUE(dec$ok)) { plot.new(); text(0.5, 0.5, dec$message); return(invisible()) }
  par_std(1, 2)
  plot(log(dec$h), log(dec$rho), pch = 20, col = COULEURS["bleu"], xlab = "log h", ylab = "log rho(h)",
       main = paste0("Échelle log-log ", titre, "\n(droite = hyperbolique)"))
  abline(lm(log(dec$rho) ~ log(dec$h)), col = COULEURS["rouge"], lwd = 2)
  plot(dec$h, log(dec$rho), pch = 20, col = COULEURS["bleu"], xlab = "h", ylab = "log rho(h)",
       main = paste0("Échelle semi-log ", titre, "\n(droite = exponentielle)"))
  abline(lm(log(dec$rho) ~ dec$h), col = COULEURS["rouge"], lwd = 2)
  par_std()
}
