# =============================================================================
# tests/test_modules.R  --  Tests de correction du code, SANS données de marché
# =============================================================================
# À lancer depuis la racine du projet :   source("tests/test_modules.R")
# Chaque test compare une brique à une référence indépendante (intégrale numérique, algèbre linéaire
# directe, identité mathématique, simulation à paramètres connus). Durée : 1 à 3 minutes.
# =============================================================================

source("R/00_utils.R", encoding = "UTF-8")
invisible(sourcer_modules("R"))
set.seed(123)

resultats <- data.frame(test = character(), valeur = character(), statut = character(), stringsAsFactors = FALSE)
verifier <- function(nom, condition, valeur = "") {
  resultats[nrow(resultats) + 1, ] <<- list(nom, as.character(valeur), ifelse(isTRUE(condition), "PASS", "FAIL"))
}

# ---- T1. Autocovariances du bruit fractionnaire FI(d) : formule récursive = intégrale spectrale -----------
d0 <- 0.3
gi <- function(h) 2 * integrate(function(l) (2 * sin(l / 2))^(-2 * d0) * cos(h * l) / (2 * pi), 0, pi,
                                 subdivisions = 2000, rel.tol = 1e-10)$value
ref <- sapply(c(0, 1, 5, 20), gi); fo <- fi_acov(d0, 20)[c(1, 2, 6, 21)]
verifier("T1  FI(d) : gamma(h) formule vs intégrale", max(abs(ref - fo)) < 1e-6, signif(max(abs(ref - fo)), 3))

# ---- T2. FARIMA(1,d,1) : convolution FI x ARMA = intégrale spectrale ------------------------------------------------
d1 <- 0.25; ar1 <- 0.5; ma1 <- 0.3
gam <- acov_farima(60, d1, ar1, ma1)
spec <- function(l) shape_spectrale(l, d1, ar1, ma1) / (2 * pi)
gint <- function(h) 2 * integrate(function(l) spec(l) * cos(h * l), 0, pi, subdivisions = 3000, rel.tol = 1e-10)$value
ref <- sapply(c(0, 1, 10, 50), gint)
verifier("T2  FARIMA(1,d,1) : acov convolution vs intégrale", max(abs(ref - gam[c(1, 2, 11, 51)])) < 1e-5,
         signif(max(abs(ref - gam[c(1, 2, 11, 51)])), 3))

# ---- T3. Vraisemblance de Durbin-Levinson = calcul gaussien direct (Cholesky) -----------------------------------------
n <- 80; gam <- acov_farima(n, d1, ar1, ma1); G <- toeplitz(gam)
x <- drop(t(chol(G)) %*% rnorm(n)) + 0.4
fin <- nll_exact_param(list(d = d1, ar = ar1, ma = ma1), x)
Gi <- solve(G); un <- rep(1, n); mu_d <- drop(un %*% Gi %*% x) / drop(un %*% Gi %*% un)
s2_d <- drop(t(x - mu_d) %*% Gi %*% (x - mu_d)) / n
val_d <- n * log(s2_d) + as.numeric(determinant(G)$modulus)
verifier("T3a Durbin-Levinson : -2 loglik concentrée vs Cholesky", abs(fin$val - val_d) < 1e-6, signif(abs(fin$val - val_d), 3))
verifier("T3b moyenne GLS profilée = GLS direct", abs(fin$mu - mu_d) < 1e-8, signif(abs(fin$mu - mu_d), 3))

# ---- T4. Convention de signes : pi * psi = identité (AR, MA et d à la fois) ---------------------------------------------
N <- 400; e <- rnorm(N)
xs <- sim_farima_innov(e, d = 0.2, ar = c(0.4, -0.2), ma = c(0.3))
e_rec <- conv_causale(xs, pi_farima(c(0.4, -0.2), 0.3, 0.2, N - 1))
verifier("T4a pi_farima inverse exactement psi_farima (filtre AR(inf) <-> MA(inf))", max(abs(e_rec - e)) < 1e-8,
         signif(max(abs(e_rec - e)), 3))
y <- rnorm(300)
verifier("T4b (1-B)^d puis (1-B)^(-d) = identité", max(abs(diff_frac(diff_frac(y, 0.35), -0.35) - y)) < 1e-9)

# ---- T5. Reparamétrisation -----------------------------------------------------------------------------------------------------
phi <- c(0.5, 0.2, -0.1); r <- ar_vers_pacf(phi)
verifier("T5a pacf <-> ar : aller-retour", max(abs(pacf_vers_ar(r) - phi)) < 1e-12)
u <- c(0.7, -1.2, 0.4, 2.0, -3.0)
pr <- u_vers_param(u, p = 2, q = 2)
rac <- verifier_racines(pr$ar, pr$ma)
verifier("T5b u -> (d, ar, ma) toujours stationnaire/inversible même pour |u| grand",
         abs(pr$d) < 0.5 && rac$stationnaire && rac$inversible)
u2 <- param_vers_u(pr$d, pr$ar, pr$ma)
verifier("T5c param <-> u : aller-retour", max(abs(u2 - u)) < 1e-6, signif(max(abs(u2 - u)), 3))

# ---- T6. Simulation Davies-Harte : autocovariances empiriques = théoriques --------------------------------------------------------
n <- 256; th <- acov_farima(n, 0.3, 0.4, -0.2)
S <- replicate(1500, sim_farima(n, 0.3, 0.4, -0.2))
emp <- sapply(c(0, 1, 5, 50), function(h) mean(S[(h + 1):n, ] * S[1:(n - h), ]))
rel <- max(abs(emp / th[c(1, 2, 6, 51)] - 1))
verifier("T6  Davies-Harte : autocov empirique vs théorique (erreur relative max < 6 %)", rel < 0.06, signif(rel, 3))

# ---- T7. Récupération de paramètres : Whittle ------------------------------------------------------------------------------------------
est <- t(replicate(25, { x <- sim_farima(2500, 0.2, 0.5, numeric(0)); coef(estim_whittle(x, 1, 0)) }))
verifier("T7a Whittle FARIMA(1,d,0) : biais sur d < 0.03", abs(mean(est[, "d"]) - 0.2) < 0.03, signif(mean(est[, "d"]), 3))
verifier("T7b Whittle : biais sur phi < 0.04", abs(mean(est[, "ar1"]) - 0.5) < 0.04, signif(mean(est[, "ar1"]), 3))
se_th <- estim_whittle(sim_farima(2500, 0.2, 0.5, numeric(0)), 1, 0)$coefs$se[1]
verifier("T7c SE asymptotique de d ~ sd empirique (rapport entre 0.7 et 1.4)",
         sd(est[, "d"]) / se_th > 0.7 && sd(est[, "d"]) / se_th < 1.4, signif(sd(est[, "d"]) / se_th, 3))
se_fi <- estim_whittle(sim_farima(3000, 0.2), 0, 0)$coefs$se
verifier("T7d FI(d) pur : se(d) = sqrt(6/(pi^2 n)) à 5 % près (discrétisation de Fourier)",
         abs(se_fi / sqrt(6 / (pi^2 * 3000)) - 1) < 0.05, signif(se_fi / sqrt(6 / (pi^2 * 3000)), 4))

# ---- T8. MLE exact et cohérence avec Whittle ---------------------------------------------------------------------------------------------
x <- sim_farima(700, 0.2, 0.5, numeric(0), mu = 0.01)
w <- estim_whittle(x, 1, 0); m <- estim_mle_exact(x, 1, 0, init = w)
verifier("T8a MLE exact vs Whittle : |d_exact - d_whittle| < 0.08", abs(m$d - w$d) < 0.08, signif(m$d - w$d, 3))
verifier("T8b MLE exact : moyenne GLS proche de 0.01 (à 4 se près)", abs(m$mu - 0.01) < 4 * se_moyenne_modele(m), signif(m$mu, 3))
verifier("T8c MLE exact : log-vraisemblance finie et >= celle d'un bruit blanc", is.finite(m$criteres["loglik"]) &&
           m$criteres["loglik"] >= estim_mle_exact(x, 0, 0, d_fixe = 0)$criteres["loglik"] - 1e-6)

# ---- T9. Prévision : AR(infini) tronqué vs prédicteur linéaire optimal exact ----------------------------------------------------------------
vp <- verifier_prevision(m, h = 10, K = 500)
verifier("T9a prévision AR(inf) vs exacte : écart max < 3 % de l'écart-type", max(abs(vp$ecart_rel_sd)) < 0.03, signif(max(abs(vp$ecart_rel_sd)), 3))
pf <- prevoir_farima(m, 10)
G500 <- toeplitz(acov_farima(510, m$d, m$ar, m$ma, m$sigma2)[1:500]); g <- acov_farima(510, m$d, m$ar, m$ma, m$sigma2)
mse_ex <- sapply(1:10, function(s) { cc <- g[(500 + s):(s + 1)]; g[1] - sum(cc * solve(G500, cc)) })
verifier("T9b MSE(h) par poids psi = MSE exacte (écart relatif < 1 %)", max(abs(pf$se^2 / mse_ex - 1)) < 0.01,
         signif(max(abs(pf$se^2 / mse_ex - 1)), 3))

# ---- T10. Estimateurs semi-paramétriques ---------------------------------------------------------------------------------------------------------
xs <- sim_farima(4000, 0.3)
lw <- est_local_whittle(xs, 0.7); gp <- est_gph(xs, 0.6)
verifier("T10a Local Whittle : d proche de 0.3 (à 3 se près)", abs(lw$d - 0.3) < 3 * lw$se, signif(lw$d, 3))
verifier("T10b GPH : d proche de 0.3 (à 3 se près)", abs(gp$d - 0.3) < 3 * gp$se, signif(gp$d, 3))
bb <- rnorm(4000)
verifier("T10c bruit blanc : Local Whittle ne rejette pas d = 0 (p > 0.01)", est_local_whittle(bb, 0.65)$p_value > 0.01)
lo <- test_lo_rs(bb)
verifier("T10d R/S de Lo : bruit blanc dans la région d'acceptation", lo$statistique > 0.809 && lo$statistique < 1.862, signif(lo$statistique, 3))

# ---- T11. Diagnostics -------------------------------------------------------------------------------------------------------------------------------------
d <- diagnostics_residus(m)
verifier("T11a résidus d'un modèle correct : Lobato (20 retards) ne rejette pas (p > 0.01)",
         d$portmanteau$p_Lobato[d$portmanteau$retards == 20] > 0.01,
         signif(d$portmanteau$p_Lobato[d$portmanteau$retards == 20], 3))
verifier("T11b ratio spectral moyen proche de 1", abs(d$spectral$valeur[1] - 1) < 0.1, signif(d$spectral$valeur[1], 3))

# ---- Bilan ---------------------------------------------------------------------------------------------------------------------------------------------------
print(resultats, row.names = FALSE, right = FALSE)
cat(sprintf("\n%d tests : %d PASS, %d FAIL\n", nrow(resultats), sum(resultats$statut == "PASS"), sum(resultats$statut == "FAIL")))
if (any(resultats$statut == "FAIL")) warning("Au moins un test a échoué : ne pas utiliser les résultats avant correction.")
invisible(resultats)
