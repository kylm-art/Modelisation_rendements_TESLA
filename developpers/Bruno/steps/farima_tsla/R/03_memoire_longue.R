# =============================================================================
# 03_memoire_longue.R  --  Détection et estimation SEMI-PARAMÉTRIQUE de la mémoire longue
# =============================================================================
# Idée : une série à mémoire longue (0 < d < 0.5) a une densité spectrale
#        f(lambda) ~ c * lambda^(-2d)  quand lambda -> 0   (pôle en 0).
# Les estimateurs semi-paramétriques n'utilisent que les m plus basses fréquences de Fourier
# lambda_j = 2*pi*j/n, j = 1..m  (m << n) et ne spécifient PAS la dynamique courte (ARMA) :
# ils sont robustes à la mauvaise spécification ARMA mais plus variables que l'estimation
# paramétrique (module 05). Ils servent de point de départ et de contrôle croisé.
#
# Estimateurs fournis :
#   - GPH (Geweke & Porter-Hudak) : régression log-périodogramme
#   - Local Whittle / GSE "Gaussian Semiparametric Estimator" (Robinson 1995)
#   - R/S de Lo (test), Hurst R/S, DFA, variance agrégée (méthodes temporelles)
# =============================================================================

# ---- Périodogramme aux fréquences de Fourier ------------------------------------------------
# Normalisation : I(lambda_j) = |sum_t x_t e^{-i lambda_j t}|^2 / n   (E[I] ~ 2*pi*f(lambda))
periodogramme <- function(x) {
  n <- length(x); x <- x - mean(x)
  m_max <- floor((n - 1) / 2)
  I <- Mod(fft(x))^2 / n
  j <- seq_len(m_max)
  list(lambda = 2 * pi * j / n, I = I[j + 1], n = n, m_max = m_max)
}

# ---- GPH ---------------------------------------------------------------------------------------------
# log I(lambda_j) = c - d * log(4 sin^2(lambda_j/2)) + u_j,   Var(u_j) = pi^2/6
est_gph <- function(x, alpha = CONFIG$alpha_bande, m = NULL) {
  per <- periodogramme(x)
  m <- min(m %||% floor(per$n^alpha), per$m_max)
  lam <- per$lambda[seq_len(m)]; I <- per$I[seq_len(m)]
  reg <- -log(4 * sin(lam / 2)^2)
  fit <- lm(log(I) ~ reg)
  d  <- unname(coef(fit)[2])
  se <- sqrt((pi^2 / 6) / sum((reg - mean(reg))^2))   # variance théorique de la régression
  data.frame(methode = "GPH", m = m, alpha = log(m) / log(per$n), d = d, se = se,
             z = d / se, p_value = 2 * pnorm(-abs(d / se)))
}

# ---- Local Whittle (GSE) -------------------------------------------------------------------------------
# Minimise  R(d) = log( (1/m) sum_j lambda_j^{2d} I_j ) - (2d/m) sum_j log lambda_j
# Loi limite : sqrt(m) (d_hat - d) -> N(0, 1/4)   ==>   se = 1 / (2 sqrt(m)).
est_local_whittle <- function(x, alpha = CONFIG$alpha_bande, m = NULL, bornes = c(-0.9, 0.9)) {
  per <- periodogramme(x)
  m <- min(m %||% floor(per$n^alpha), per$m_max)
  lam <- per$lambda[seq_len(m)]; I <- per$I[seq_len(m)]
  slog <- mean(log(lam))
  obj <- function(d) log(mean(lam^(2 * d) * I)) - 2 * d * slog
  d <- optimize(obj, bornes)$minimum
  se <- 1 / (2 * sqrt(m))
  data.frame(methode = "Local Whittle (GSE)", m = m, alpha = log(m) / log(per$n), d = d, se = se,
             z = d / se, p_value = 2 * pnorm(-abs(d / se)))
}

# ---- Sensibilité à la bande passante m = n^alpha ------------------------------------------------------------
# Une estimation fiable de d doit être STABLE quand alpha varie ; une dérive systématique trahit
# soit un biais de la dynamique courte (alpha trop grand), soit de la variance (alpha trop petit),
# soit une fausse mémoire longue (ruptures, non-stationnarité).
scan_bande_passante <- function(x, alphas = seq(0.40, 0.80, by = 0.05)) {
  do.call(rbind, lapply(alphas, function(a) {
    g <- est_gph(x, alpha = a); l <- est_local_whittle(x, alpha = a)
    data.frame(alpha = a, m = g$m, gph_d = g$d, gph_se = g$se, lw_d = l$d, lw_se = l$se)
  }))
}

plot_scan_d <- function(sc, titre = "") {
  rng <- range(c(sc$gph_d - 2 * sc$gph_se, sc$gph_d + 2 * sc$gph_se,
                 sc$lw_d - 2 * sc$lw_se, sc$lw_d + 2 * sc$lw_se, 0))
  par_std()
  plot(sc$alpha, sc$gph_d, type = "b", pch = 19, col = COULEURS["bleu"], ylim = rng,
       xlab = "alpha  (m = n^alpha)", ylab = "d estimé",
       main = paste("Stabilité de d selon la bande passante", titre))
  segments(sc$alpha, sc$gph_d - 2 * sc$gph_se, sc$alpha, sc$gph_d + 2 * sc$gph_se, col = COULEURS["bleu"])
  lines(sc$alpha + 0.004, sc$lw_d, type = "b", pch = 17, col = COULEURS["rouge"])
  segments(sc$alpha + 0.004, sc$lw_d - 2 * sc$lw_se, sc$alpha + 0.004, sc$lw_d + 2 * sc$lw_se, col = COULEURS["rouge"])
  abline(h = 0, lty = 2, col = "grey40")
  legend("topright", c("GPH (±2 se)", "Local Whittle (±2 se)"), col = COULEURS[c("bleu", "rouge")],
         pch = c(19, 17), bty = "n", cex = 0.85)
}

# ---- Log-périodogramme avec droite de pente -2d --------------------------------------------------------------
plot_periodogramme <- function(x, alpha = CONFIG$alpha_bande, titre = "") {
  per <- periodogramme(x); m <- floor(per$n^alpha)
  gph <- est_gph(x, alpha = alpha)
  reg <- -log(4 * sin(per$lambda / 2)^2)
  fit <- lm(log(per$I[seq_len(m)]) ~ reg[seq_len(m)])
  par_std()
  plot(log(per$lambda), log(per$I), pch = 20, cex = 0.35, col = "grey65",
       xlab = "log(lambda)", ylab = "log I(lambda)", main = paste("Log-périodogramme", titre))
  points(log(per$lambda[seq_len(m)]), log(per$I[seq_len(m)]), pch = 20, cex = 0.6, col = COULEURS["bleu"])
  # droite GPH : log I = c - d * log(4 sin^2(lambda/2))  ~  c + 2 d * (-log lambda) ... pente en log(lambda) = -2d
  lines(log(per$lambda[seq_len(m)]), coef(fit)[1] + coef(fit)[2] * reg[seq_len(m)], col = COULEURS["rouge"], lwd = 2)
  legend("topright", c("fréquences utilisées (m)", paste0("pente -2d, d = ", round(gph$d, 3))),
         col = COULEURS[c("bleu", "rouge")], pch = c(20, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n", cex = 0.85)
}

# ---- Test R/S modifié de Lo (1991) --------------------------------------------------------------------------------
# H0 : mémoire courte.  V = (max_k S_k - min_k S_k) / (sigma_q sqrt(n)), avec sigma_q^2 estimateur
# de variance de long terme (noyau de Bartlett, q retards) qui corrige la dépendance de courte portée.
# Région d'acceptation à 5 % : V dans [0.809 ; 1.862]  (1 % : [0.721 ; 2.098]).
test_lo_rs <- function(x, q = NULL) {
  n <- length(x); xc <- x - mean(x)
  q <- q %||% floor(4 * (n / 100)^0.25)
  S <- cumsum(xc)
  R <- max(c(0, S)) - min(c(0, S))
  gam <- function(j) sum(xc[(j + 1):n] * xc[1:(n - j)]) / n
  s2 <- gam(0)
  if (q > 0) s2 <- s2 + 2 * sum(vapply(seq_len(q), function(j) (1 - j / (q + 1)) * gam(j), numeric(1)))
  V <- R / sqrt(n * s2)
  data.frame(test = paste0("R/S modifié de Lo (q = ", q, ")"), statistique = V,
             borne_inf_5pc = 0.809, borne_sup_5pc = 1.862,
             conclusion = ifelse(V < 0.809 | V > 1.862, "mémoire courte rejetée (5 %)", "mémoire courte non rejetée"))
}

# ---- Estimateurs temporels de l'exposant de Hurst (H = d + 1/2) ------------------------------------------------------------
tailles_log <- function(n, min_s, max_s, k = 15) unique(floor(exp(seq(log(min_s), log(max_s), length.out = k))))

hurst_rs <- function(x, tailles = NULL) {
  n <- length(x); tailles <- tailles %||% tailles_log(n, 16, n / 2)
  rs <- vapply(tailles, function(s) {
    k <- floor(n / s)
    mean(vapply(seq_len(k), function(i) {
      seg <- x[((i - 1) * s + 1):(i * s)]; y <- cumsum(seg - mean(seg))
      sdv <- sqrt(mean((seg - mean(seg))^2))
      if (sdv > 0) (max(y) - min(y)) / sdv else NA_real_
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))
  fit <- lm(log(rs) ~ log(tailles)); H <- unname(coef(fit)[2])
  list(H = H, d = H - 0.5, tailles = tailles, stat = rs, fit = fit, nom = "Hurst R/S")
}

hurst_dfa <- function(x, tailles = NULL, ordre = 1) {
  y <- cumsum(x - mean(x)); n <- length(y)
  tailles <- tailles %||% tailles_log(n, 10, n / 5)
  Fs <- vapply(tailles, function(s) {
    k <- floor(n / s)
    X <- outer(seq_len(s) / s, 0:ordre, "^"); Q <- qr.Q(qr(X))
    f2 <- vapply(seq_len(k), function(i) {
      seg <- y[((i - 1) * s + 1):(i * s)]
      res <- seg - drop(Q %*% crossprod(Q, seg)); mean(res^2)
    }, numeric(1))
    sqrt(mean(f2))
  }, numeric(1))
  fit <- lm(log(Fs) ~ log(tailles)); a <- unname(coef(fit)[2])   # pour un bruit stationnaire : a = H
  list(H = a, d = a - 0.5, tailles = tailles, stat = Fs, fit = fit, nom = "DFA")
}

hurst_varagg <- function(x, tailles = NULL) {
  n <- length(x); tailles <- tailles %||% tailles_log(n, 4, n / 10)
  v <- vapply(tailles, function(m) {
    k <- floor(n / m); var(colMeans(matrix(x[seq_len(k * m)], nrow = m)))
  }, numeric(1))
  fit <- lm(log(v) ~ log(tailles)); H <- 1 + unname(coef(fit)[2]) / 2   # Var(X^(m)) ~ m^(2H-2)
  list(H = H, d = H - 0.5, tailles = tailles, stat = v, fit = fit, nom = "Variance agrégée")
}

plot_hurst <- function(rs, dfa, va) {
  par_std(1, 3)
  for (h in list(rs, dfa, va)) {
    plot(log(h$tailles), log(h$stat), pch = 19, col = COULEURS["bleu"],
         xlab = "log(taille)", ylab = "log(statistique)",
         main = paste0(h$nom, " : H = ", round(h$H, 3), ", d = ", round(h$d, 3)))
    abline(h$fit, col = COULEURS["rouge"], lwd = 2)
  }
  par_std()
}

# ---- Tableau de synthèse : toutes les estimations de d -----------------------------------------------------------------------
resume_memoire_longue <- function(x, alpha = CONFIG$alpha_bande) {
  gph1 <- est_gph(x, alpha = 0.50); gph2 <- est_gph(x, alpha = alpha)
  lw1  <- est_local_whittle(x, alpha = 0.50); lw2 <- est_local_whittle(x, alpha = alpha)
  rs <- hurst_rs(x); dfa <- hurst_dfa(x); va <- hurst_varagg(x)
  out <- rbind(
    transform(gph1, methode = "GPH (alpha = 0.50)"),
    transform(gph2, methode = paste0("GPH (alpha = ", alpha, ")")),
    transform(lw1,  methode = "Local Whittle (alpha = 0.50)"),
    transform(lw2,  methode = paste0("Local Whittle (alpha = ", alpha, ")")))[, c("methode", "d", "se", "z", "p_value")]
  tmp <- data.frame(methode = c("Hurst R/S (H - 1/2)", "DFA-1 (H - 1/2)", "Variance agrégée (H - 1/2)"),
                    d = c(rs$d, dfa$d, va$d), se = NA_real_, z = NA_real_, p_value = NA_real_)
  rbind(out, tmp)
}

# ---- Mémoire longue des transformations (r, |r|, r^2) ------------------------------------------------------------------------------
memoire_transformations <- function(y, alpha = CONFIG$alpha_bande) {
  f <- function(x, nom) { l <- est_local_whittle(x, alpha); g <- est_gph(x, alpha)
    data.frame(serie = nom, d_GPH = g$d, se_GPH = g$se, d_LW = l$d, se_LW = l$se,
               p_LW_d0 = l$p_value) }
  rbind(f(y, "r"), f(abs(y), "|r|"), f(y^2, "r^2"), f(log(abs(y) + 1e-6), "log(|r|)"))
}

# ---- Stabilité de d dans le temps (détection de FAUSSE mémoire longue) ----------------------------------------------------------------------
# Des ruptures de moyenne ou de variance peuvent imiter une mémoire longue. Si d est stable
# sur des sous-périodes, la mémoire longue est une propriété de la dynamique ; s'il varie fortement
# ou s'effondre hors des périodes de rupture, méfiance.
d_glissant <- function(x, dates = NULL, largeur = 1000, pas = 250, alpha = 0.65) {
  n <- length(x); debuts <- seq(1, n - largeur + 1, by = pas)
  res <- do.call(rbind, lapply(debuts, function(i) {
    seg <- x[i:(i + largeur - 1)]; l <- est_local_whittle(seg, alpha)
    data.frame(debut = i, fin = i + largeur - 1, d = l$d, se = l$se)
  }))
  if (!is.null(dates)) { res$date_debut <- dates[res$debut]; res$date_fin <- dates[res$fin] }
  res
}

plot_d_glissant <- function(res, titre = "") {
  mid <- (res$debut + res$fin) / 2
  par_std()
  plot(mid, res$d, type = "b", pch = 19, col = COULEURS["bleu"],
       ylim = range(c(res$d - 2 * res$se, res$d + 2 * res$se, 0)),
       xlab = "milieu de fenêtre (indice)", ylab = "d (local Whittle)", main = paste("d sur fenêtres glissantes", titre))
  segments(mid, res$d - 2 * res$se, mid, res$d + 2 * res$se, col = COULEURS["bleu"])
  abline(h = 0, lty = 2, col = "grey40")
}
