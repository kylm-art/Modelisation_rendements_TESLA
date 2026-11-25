# =============================================================================
# 06_selection.R  --  Sélection de l'ordre (p, q) sur une grille, FARIMA vs ARMA
# =============================================================================
# Principe : on estime tous les FARIMA(p, d, q) avec p <= P, q <= Q, puis on compare des critères
# d'information. Le BIC est consistant (sélectionne le vrai ordre quand n -> inf) ; l'AIC est
# plus généreux ; HQ est intermédiaire. On regarde AUSSI la proximité des racines AR/MA : deux racines
# quasi identiques signifient qu'un facteur commun se simplifie (modèle sur-paramétré, SE explosifs).
# =============================================================================

grille_farima <- function(x, P = CONFIG$P_max, Q = CONFIG$Q_max, methode = c("whittle", "exact"),
                          d_fixe = NULL, progression = FALSE) {
  methode <- match.arg(methode)
  lignes <- list(); fits <- list()
  for (p in 0:P) for (q in 0:Q) {
    if (progression) message(sprintf("  FARIMA(%d,%s,%d) [%s]", p, ifelse(is.null(d_fixe), "d", d_fixe), q, methode))
    f <- tryCatch(estim_farima(x, p, q, methode, d_fixe), error = function(e) NULL)
    if (is.null(f)) next
    rac <- verifier_racines(f$ar, f$ma)
    nm <- sprintf("%d,%d", p, q); fits[[nm]] <- f
    lignes[[nm]] <- data.frame(p = p, q = q, d = f$d, sigma2 = unname(f$criteres["sigma2"]),
                               AIC = unname(f$criteres["AIC"]), BIC = unname(f$criteres["BIC"]),
                               HQ = unname(f$criteres["HQ"]), k = unname(f$criteres["k"]),
                               converge = f$converge, dist_racines = rac$distance_min_ar_ma,
                               stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, lignes); rownames(tab) <- NULL
  tab$racines_quasi_communes <- !is.na(tab$dist_racines) & tab$dist_racines < 0.10
  list(table = tab, fits = fits, methode = methode, d_fixe = d_fixe)
}

# Classement par critère avec écart au meilleur (delta). Règle usuelle : delta < 2 = équivalent.
classer_modeles <- function(res, critere = c("BIC", "AIC", "HQ")) {
  critere <- match.arg(critere)
  tab <- res$table[order(res$table[[critere]]), ]
  tab$delta <- tab[[critere]] - min(tab[[critere]])
  rownames(tab) <- NULL
  tab
}

meilleur_modele <- function(res, critere = c("BIC", "AIC", "HQ"), exclure_redondants = TRUE) {
  critere <- match.arg(critere)
  tab <- classer_modeles(res, critere)
  if (exclure_redondants) tab <- tab[!tab$racines_quasi_communes, , drop = FALSE]
  tab <- tab[tab$converge, , drop = FALSE]
  res$fits[[sprintf("%d,%d", tab$p[1], tab$q[1])]]
}

plot_grille_critere <- function(res, critere = "BIC") {
  tab <- res$table; P <- max(tab$p); Q <- max(tab$q)
  M <- matrix(NA_real_, P + 1, Q + 1)
  for (i in seq_len(nrow(tab))) M[tab$p[i] + 1, tab$q[i] + 1] <- tab[[critere]][i]
  delta <- M - min(M, na.rm = TRUE)
  par_std()
  image(0:Q, 0:P, t(delta), col = hcl.colors(30, "Blues 3", rev = FALSE), xlab = "q (MA)", ylab = "p (AR)",
        axes = FALSE, main = paste0(critere, " - min  (", res$methode, ")"))
  axis(1, 0:Q); axis(2, 0:P)
  for (p in 0:P) for (q in 0:Q) if (!is.na(delta[p + 1, q + 1]))
    text(q, p, round(delta[p + 1, q + 1], 1), cex = 0.9,
         font = ifelse(delta[p + 1, q + 1] == 0, 2, 1))
}

# FARIMA contre ARMA : on refait la grille avec d = 0 imposé et on compare les meilleurs de chaque famille.
comparer_farima_arma <- function(x, P = CONFIG$P_max, Q = CONFIG$Q_max, methode = "whittle", critere = "BIC",
                                 g1 = NULL) {
  g1 <- g1 %||% grille_farima(x, P, Q, methode)      # g1 : grille FARIMA déjà calculée (évite de la refaire)
  g0 <- grille_farima(x, P, Q, methode, d_fixe = 0)
  b1 <- classer_modeles(g1, critere)[1, ]; b0 <- classer_modeles(g0, critere)[1, ]
  out <- rbind(data.frame(famille = "FARIMA (d libre)", p = b1$p, q = b1$q, d = b1$d, critere = b1[[critere]]),
               data.frame(famille = "ARMA (d = 0)",     p = b0$p, q = b0$q, d = 0,     critere = b0[[critere]]))
  names(out)[5] <- critere
  list(table = out, grille_farima = g1, grille_arma = g0)
}
