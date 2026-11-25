# =============================================================================
# 09_backtest.R  --  Évaluation hors échantillon (pseudo temps réel) et comparaison de modèles
# =============================================================================
# Protocole :  pour chaque origine t (de n_init à n-1), on n'utilise QUE x_1..x_t, on prévoit x_{t+h},
# puis on compare à la valeur réalisée. Les paramètres sont ré-estimés tous les `pas_refit` jours.
# Aucune information du futur n'entre dans une prévision (pas de fuite).
#
# Un "modèle" est une liste de deux fonctions :
#     $ajuster(train)           -> objet (paramètres estimés)
#     $prevoir(obj, train, h)   -> vecteur de longueur h
# ce qui permet d'ajouter n'importe quel concurrent sans toucher au reste du code.
# =============================================================================

# ---- Fabrique de modèles --------------------------------------------------------------------------------------------
modele_farima <- function(p, q, d_fixe = NULL, methode = "whittle", nom = NULL, K = 1000) {
  nom <- nom %||% sprintf("%s(%d,%s,%d)", ifelse(is.null(d_fixe), "FARIMA", "ARMA"), p,
                          ifelse(is.null(d_fixe), "d", d_fixe), q)
  list(nom = nom,
       ajuster = function(train) estim_farima(train, p, q, methode, d_fixe),
       prevoir = function(obj, train, h) {
         obj$x <- train; obj$n <- length(train); obj$mu <- mean(train)     # paramètres figés, données mises à jour
         prevoir_farima(obj, h, niveaux = numeric(0), K = min(K, length(train)))$prevision
       })
}
modele_moyenne <- function(nom = "Moyenne historique")
  list(nom = nom, ajuster = function(train) NULL, prevoir = function(obj, train, h) rep(mean(train), h))
modele_zero <- function(nom = "Zéro (marche aléatoire)")
  list(nom = nom, ajuster = function(train) NULL, prevoir = function(obj, train, h) rep(0, h))
modele_mm <- function(k = 22, nom = NULL)
  list(nom = nom %||% paste0("Moyenne mobile ", k, " j"), ajuster = function(train) NULL,
       prevoir = function(obj, train, h) rep(mean(tail(train, k)), h))
modele_ewma <- function(lambda = 0.94, nom = NULL)
  list(nom = nom %||% paste0("EWMA (", lambda, ")"), ajuster = function(train) NULL,
       prevoir = function(obj, train, h) {
         niv <- as.numeric(stats::filter((1 - lambda) * train, lambda, method = "recursive"))
         rep(niv[length(niv)], h) })

# ---- Moteur de backtest ---------------------------------------------------------------------------------------------------
backtest <- function(x, modeles, h_set = c(1, 5), n_init = floor(0.6 * length(x)), pas_refit = 20,
                     fenetre = c("expansion", "glissante"), largeur = NULL, progression = TRUE) {
  fenetre <- match.arg(fenetre)
  if (fenetre == "glissante" && is.null(largeur)) stop("fenetre = 'glissante' exige `largeur`.", call. = FALSE)
  n <- length(x); H <- max(h_set); origines <- n_init:(n - 1)
  noms <- vapply(modeles, function(m) m$nom, character(1)); names(modeles) <- noms
  prev <- lapply(noms, function(nm) matrix(NA_real_, length(origines), length(h_set),
                                           dimnames = list(NULL, paste0("h", h_set))))
  names(prev) <- noms
  objets <- vector("list", length(modeles)); names(objets) <- noms

  for (i in seq_along(origines)) {
    t <- origines[i]
    train <- x[(if (fenetre == "glissante") max(1, t - largeur + 1) else 1):t]
    if (i == 1 || (i - 1) %% pas_refit == 0) {
      if (progression && (i - 1) %% (5 * pas_refit) == 0) message(sprintf("  backtest : origine %d / %d", i, length(origines)))
      for (nm in noms) objets[[nm]] <- tryCatch(modeles[[nm]]$ajuster(train), error = function(e) NULL)
    }
    for (nm in noms) {
      f <- tryCatch(modeles[[nm]]$prevoir(objets[[nm]], train, H), error = function(e) rep(NA_real_, H))
      prev[[nm]][i, ] <- f[h_set]
    }
  }
  reel <- vapply(h_set, function(h) ifelse(origines + h <= n, x[pmin(origines + h, n)], NA_real_),
                 numeric(length(origines)))
  colnames(reel) <- paste0("h", h_set)
  list(origines = origines, reel = reel, prev = prev, h_set = h_set, noms = noms, n = n)
}

# ---- Tests de précision prédictive ----------------------------------------------------------------------------------------------
# Variance de long terme tronquée à h-1 retards (Diebold-Mariano)
vlt_dm <- function(d, h) {
  n <- length(d); dc <- d - mean(d)
  g <- vapply(0:(h - 1), function(k) sum(dc[(k + 1):n] * dc[1:(n - k)]) / n, numeric(1))
  g[1] + 2 * sum(g[-1])
}

# Diebold-Mariano avec correction de Harvey-Leybourne-Newbold. d_t = L(e_ref) - L(e_mod) ; H1 : le modèle fait mieux.
dm_test <- function(e_ref, e_mod, h = 1, puissance = 2) {
  d <- abs(e_ref)^puissance - abs(e_mod)^puissance; n <- length(d)
  V <- vlt_dm(d, h); if (!is.finite(V) || V <= 0) return(c(stat = NA_real_, p_value = NA_real_))
  dm <- mean(d) / sqrt(V / n)
  hln <- dm * sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  c(stat = hln, p_value = pt(hln, n - 1, lower.tail = FALSE))
}

# Clark-West (2007) : version adaptée aux modèles EMBOÎTÉS (ex. FARIMA contre moyenne constante),
# où le DM est conservateur. f_t = e_ref^2 - [e_mod^2 - (f_ref - f_mod)^2] ; H1 : le modèle fait mieux.
cw_test <- function(y, f_ref, f_mod, h = 1) {
  ft <- (y - f_ref)^2 - ((y - f_mod)^2 - (f_ref - f_mod)^2); n <- length(ft)
  V <- vlt_dm(ft, h); if (!is.finite(V) || V <= 0) return(c(stat = NA_real_, p_value = NA_real_))
  s <- mean(ft) / sqrt(V / n)
  c(stat = s, p_value = pnorm(s, lower.tail = FALSE))
}

evaluer_backtest <- function(bt, ref = bt$noms[1], h = 1) {
  col <- match(h, bt$h_set); if (is.na(col)) stop("h absent du backtest.", call. = FALSE)
  y <- bt$reel[, col]
  ok <- !is.na(y) & Reduce(`&`, lapply(bt$prev, function(m) !is.na(m[, col])))
  y <- y[ok]; P <- lapply(bt$prev, function(m) m[ok, col]); E <- lapply(P, function(f) y - f)
  mse_ref <- mean(E[[ref]]^2); base <- max(mean(y > 0), mean(y < 0))
  tab <- do.call(rbind, lapply(bt$noms, function(nm) {
    e <- E[[nm]]; f <- P[[nm]]
    dm <- if (nm == ref) c(stat = NA, p_value = NA) else dm_test(E[[ref]], e, h)
    cw <- if (nm == ref) c(stat = NA, p_value = NA) else cw_test(y, P[[ref]], f, h)
    nz <- y != 0 & f != 0; hit <- if (sum(nz) > 0) mean(sign(f[nz]) == sign(y[nz])) else NA_real_
    p_hit <- if (sum(nz) > 0) binom.test(sum(sign(f[nz]) == sign(y[nz])), sum(nz), p = base, alternative = "greater")$p.value else NA_real_
    data.frame(modele = nm, n_prev = length(y), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)), biais = mean(e),
               R2_OOS = 1 - mean(e^2) / mse_ref, DM_stat = dm["stat"], DM_p = dm["p_value"],
               CW_stat = cw["stat"], CW_p = cw["p_value"], hit_rate = hit, hit_p = p_hit, row.names = NULL)
  }))
  tab
}

evaluer_tous_horizons <- function(bt, ref = bt$noms[1]) {
  do.call(rbind, lapply(bt$h_set, function(h) cbind(horizon = h, evaluer_backtest(bt, ref, h))))
}

# ---- Graphique : différence cumulée des erreurs quadratiques (CSFE) --------------------------------------------------------------------
# Pente positive = le modèle bat la référence à cet instant ; une courbe qui monte puis redescend
# révèle une supériorité instable.
plot_csfe <- function(bt, ref = bt$noms[1], h = 1) {
  col <- match(h, bt$h_set); y <- bt$reel[, col]; ok <- !is.na(y)
  e_ref <- (y - bt$prev[[ref]][, col])[ok]
  autres <- setdiff(bt$noms, ref); cols <- unname(COULEURS)[seq_along(autres)]
  cs <- lapply(autres, function(nm) cumsum(e_ref^2 - ((y - bt$prev[[nm]][, col])[ok])^2))
  par_std()
  plot(bt$origines[ok], cs[[1]], type = "n", ylim = range(unlist(cs), na.rm = TRUE), xlab = "origine (indice)",
       ylab = "somme cumulée  e_ref^2 - e_mod^2", main = paste0("Gain cumulé contre ", ref, "  (h = ", h, ")"))
  for (i in seq_along(autres)) lines(bt$origines[ok], cs[[i]], col = cols[i], lwd = 2)
  abline(h = 0, lty = 2)
  legend("topleft", autres, col = cols, lwd = 2, bty = "n", cex = 0.8)
}
