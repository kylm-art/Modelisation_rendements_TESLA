# =============================================================================
# 00_utils.R  --  Configuration globale, chargement des packages, utilitaires
# =============================================================================
# Convention du projet (à garder en tête dans tous les modules) :
#
#   Le modèle FARIMA(p, d, q) s'écrit
#
#        Phi(B) (1 - B)^d (X_t - mu) = Theta(B) eps_t
#
#   avec  Phi(B)   = 1 - phi_1 B - ... - phi_p B^p      (partie AR)
#         Theta(B) = 1 + theta_1 B + ... + theta_q B^q  (partie MA)
#   et   eps_t bruit blanc de variance sigma2.
#   C'est la convention de stats::arima : ar = (phi_1..phi_p), ma = (theta_1..theta_q).
#
#   "Gauss" désigne ici l'estimation par quasi-vraisemblance gaussienne :
#     - Whittle (vraisemblance gaussienne dans le domaine fréquentiel),
#     - maximum de vraisemblance gaussien EXACT (Durbin-Levinson),
#     - et, en semi-paramétrique, l'estimateur gaussien de Robinson (local Whittle / GSE).
# =============================================================================

# ---- Configuration -----------------------------------------------------------
CONFIG <- list(
  ticker        = "TSLA",
  debut         = "2015-01-01",
  fin           = "2026-11-23",
  fichier_cache = "data/TSLA.csv",   # copie locale (écrite après un téléchargement réussi)
  graine        = 20261125,
  jours_an      = 252,
  P_max         = 2,                 # ordre AR maximal exploré dans la grille
  Q_max         = 2,                 # ordre MA maximal exploré dans la grille
  alpha_bande   = 0.65,              # m = n^alpha pour les estimateurs semi-paramétriques
  d_max         = 0.49               # borne de stationnarité imposée à d (|d| < 0.5)
)

# ---- Chargement des packages --------------------------------------------------
# requis    : le projet ne tourne pas sans eux
# optionnels: étapes de contrôle croisé / extensions (le code se dégrade proprement)
PACKAGES_REQUIS      <- c("quantmod", "xts", "zoo", "tseries")
PACKAGES_OPTIONNELS  <- c("fracdiff", "rugarch")

charger_packages <- function(requis = PACKAGES_REQUIS, optionnels = PACKAGES_OPTIONNELS) {
  manquants <- requis[!vapply(requis, requireNamespace, logical(1), quietly = TRUE)]
  if (length(manquants)) {
    stop("Packages requis manquants : ", paste(manquants, collapse = ", "),
         "\n  -> install.packages(c(", paste0('"', manquants, '"', collapse = ", "), "))",
         call. = FALSE)
  }
  suppressPackageStartupMessages(
    for (p in requis) library(p, character.only = TRUE)
  )
  dispo <- vapply(optionnels, requireNamespace, logical(1), quietly = TRUE)
  names(dispo) <- optionnels
  invisible(dispo)
}

pkg_dispo <- function(p) requireNamespace(p, quietly = TRUE)

# ---- Petits utilitaires --------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

fmt_p <- function(p, digits = 4) {
  ifelse(is.na(p), "NA",
         ifelse(p < 10^(-digits), paste0("< ", format(10^(-digits), scientific = FALSE)),
                formatC(p, digits = digits, format = "f")))
}

arrondir_df <- function(df, digits = 4) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  df
}

# Verdict textuel (affiché dans le rapport) à partir d'une p-value
verdict_test <- function(p, alpha = 0.05, rejet = "H0 rejetée", non_rejet = "H0 non rejetée") {
  if (is.na(p)) return("indéterminé")
  if (p < alpha) rejet else non_rejet
}

# Interprétation automatique d'un paramètre de mémoire longue d
interpreter_d <- function(d, se = NA_real_, niveau = 0.05) {
  z <- if (is.na(se)) NA_real_ else d / se
  signif <- !is.na(z) && abs(z) > qnorm(1 - niveau / 2)
  regime <-
    if (d >= 0.5)      "non stationnaire (d >= 0,5)"
    else if (d > 0.05) "mémoire longue stationnaire (0 < d < 0,5) : autocorrélations hyperboliques, non sommables"
    else if (d < -0.05) "antipersistance (d < 0) : densité spectrale nulle en 0 (sur-différenciation possible)"
    else               "mémoire courte (d proche de 0) : autocorrélations à décroissance exponentielle"
  paste0("d = ", formatC(d, digits = 3, format = "f"),
         if (!is.na(se)) paste0(" (se = ", formatC(se, digits = 3, format = "f"), ")") else "",
         " -> ", regime,
         if (!is.na(z)) paste0(" ; ", if (signif) "significativement" else "non significativement",
                               " différent de 0 au seuil de ", niveau * 100, " %") else "")
}

# ---- Graphiques : réglages cohérents -------------------------------------------
par_std <- function(nr = 1, nc = 1) {
  par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), mgp = c(2.4, 0.7, 0), cex.main = 1.0)
}

COULEURS <- c(bleu = "#1f77b4", rouge = "#d62728", vert = "#2ca02c",
              orange = "#ff7f0e", gris = "grey55", violet = "#9467bd")

# ---- Sources des modules --------------------------------------------------------
# Les fichiers sont sourcés dans l'ordre numérique (chaque module n'utilise que les précédents).
sourcer_modules <- function(dossier = "R") {
  fichiers <- sort(list.files(dossier, pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
  fichiers <- fichiers[basename(fichiers) != "00_utils.R"]
  for (f in fichiers) source(f, encoding = "UTF-8")
  invisible(basename(fichiers))
}
