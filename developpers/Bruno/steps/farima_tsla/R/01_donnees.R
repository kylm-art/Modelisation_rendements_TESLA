# =============================================================================
# 01_donnees.R  --  Acquisition, nettoyage et contrôle qualité des données
# =============================================================================
# Sortie principale : une liste "jeu" contenant
#   $ohlc   : xts des prix bruts (Open, High, Low, Close, Volume, Adjusted)
#   $prix   : xts du prix ajusté (splits + dividendes)
#   $r      : xts des rendements logarithmiques  r_t = log(P_t / P_{t-1})
#   $y      : vecteur numérique des rendements (entrée des modules suivants)
#   $dates  : dates associées à $y
#   $source : "yahoo" ou "cache"
# =============================================================================

# ---- 1. Téléchargement avec repli sur un fichier local -------------------------
charger_donnees <- function(ticker = CONFIG$ticker, debut = CONFIG$debut, fin = CONFIG$fin,
                            fichier_cache = CONFIG$fichier_cache, forcer_cache = FALSE) {
  ohlc <- NULL; source <- NA_character_

  if (!forcer_cache) {
    ohlc <- tryCatch(
      suppressWarnings(getSymbols(ticker, src = "yahoo", from = debut, to = fin,
                                  auto.assign = FALSE)),
      error = function(e) { message("Téléchargement impossible : ", conditionMessage(e)); NULL })
    if (!is.null(ohlc)) {
      source <- "yahoo"
      dir.create(dirname(fichier_cache), showWarnings = FALSE, recursive = TRUE)
      zoo::write.zoo(ohlc, fichier_cache, sep = ",")
    }
  }
  if (is.null(ohlc)) {
    if (!file.exists(fichier_cache))
      stop("Aucune donnée : pas de réseau et pas de fichier '", fichier_cache, "'.", call. = FALSE)
    z <- zoo::read.zoo(fichier_cache, header = TRUE, sep = ",", format = "%Y-%m-%d")
    ohlc <- xts::as.xts(z)
    source <- "cache"
  }
  attr(ohlc, "source") <- source
  ohlc
}

# ---- 2. Construction des rendements ---------------------------------------------
construire_jeu <- function(ohlc) {
  prix <- quantmod::Ad(ohlc)
  colnames(prix) <- "prix"
  r <- na.omit(diff(log(prix)))
  colnames(r) <- "r"
  list(ohlc = ohlc, prix = prix, r = r,
       y = as.numeric(r), dates = zoo::index(r),
       source = attr(ohlc, "source"))
}

# ---- 3. Contrôle qualité ----------------------------------------------------------
# Renvoie un data.frame (test, résultat, statut) + détails. Aucune correction automatique :
# on SIGNALE, on ne modifie pas silencieusement les données.
controle_qualite <- function(jeu, seuil_extreme = 5, trou_max_jours = 5) {
  p <- as.numeric(jeu$prix); d <- jeu$dates; y <- jeu$y
  z <- (y - mean(y)) / sd(y)
  ecarts <- as.numeric(diff(as.Date(zoo::index(jeu$prix))))
  extremes <- which(abs(z) > seuil_extreme)

  res <- data.frame(
    controle = c("Nombre d'observations (rendements)",
                 "Période couverte",
                 "Valeurs manquantes dans les prix",
                 "Prix négatifs ou nuls",
                 "Dates dupliquées",
                 "Dates non croissantes",
                 paste0("Trous > ", trou_max_jours, " jours calendaires"),
                 "Rendements exactement nuls",
                 paste0("Rendements extrêmes (|z| > ", seuil_extreme, ")"),
                 "Rendement minimal / maximal"),
    resultat = c(length(y),
                 paste(range(as.character(d)), collapse = " -> "),
                 sum(is.na(p)),
                 sum(p <= 0, na.rm = TRUE),
                 sum(duplicated(zoo::index(jeu$prix))),
                 sum(ecarts <= 0),
                 sum(ecarts > trou_max_jours),
                 sum(y == 0),
                 length(extremes),
                 paste(round(range(y), 4), collapse = " / ")),
    stringsAsFactors = FALSE)
  res$statut <- c("info", "info",
                  ifelse(sum(is.na(p)) == 0, "OK", "ALERTE"),
                  ifelse(sum(p <= 0, na.rm = TRUE) == 0, "OK", "ALERTE"),
                  ifelse(sum(duplicated(zoo::index(jeu$prix))) == 0, "OK", "ALERTE"),
                  ifelse(sum(ecarts <= 0) == 0, "OK", "ALERTE"),
                  ifelse(sum(ecarts > trou_max_jours) == 0, "OK", "à examiner"),
                  ifelse(sum(y == 0) / length(y) < 0.02, "OK", "à examiner"),
                  ifelse(length(extremes) == 0, "OK", "à examiner"),
                  "info")

  list(table = res,
       extremes = data.frame(date = d[extremes], rendement = y[extremes], z = z[extremes]),
       trous = data.frame(date = zoo::index(jeu$prix)[-1][ecarts > trou_max_jours],
                          jours = ecarts[ecarts > trou_max_jours]))
}

# ---- 4. Découpage apprentissage / test (sans fuite d'information) --------------------
decouper_echantillon <- function(y, part_test = 0.2) {
  n <- length(y); n_test <- floor(n * part_test)
  list(app = y[seq_len(n - n_test)], test = y[(n - n_test + 1):n],
       n_app = n - n_test, n_test = n_test)
}
