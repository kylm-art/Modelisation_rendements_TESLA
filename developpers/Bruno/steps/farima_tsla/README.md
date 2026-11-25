# Pipeline FARIMA gaussien — rendements TESLA

Modélisation FARIMA(p, d, q) par quasi-vraisemblance gaussienne (Whittle, MLE exact, local Whittle/GSE),
de l'exploration des données aux prévisions et à leur validation hors échantillon. Code modulaire en R,
rapport Quarto commenté.

## Lancer

```r
# 1) packages
install.packages(c("quantmod", "xts", "zoo", "tseries"))   # requis
install.packages(c("fracdiff", "rugarch"))                  # optionnels (contrôle croisé, FARIMA-GARCH)

# 2) vérifier le code (1 à 3 min, aucune donnée de marché nécessaire)
source("tests/test_modules.R")

# 3) produire le rapport
quarto::quarto_render("analyse_FARIMA_TSLA.qmd")     # ou : quarto render analyse_FARIMA_TSLA.qmd
```

Dans le `.qmd`, `MODE_RAPIDE <- TRUE` réduit la grille, le bootstrap et le backtest pour un premier passage
(quelques minutes) ; `EXECUTER_TESTS <- TRUE` intègre les tests au rapport. Le rendu complet est long
(MLE exact en O(n²), bootstrap, backtests) : les chunks lourds sont mis en cache (`cache: true`).

Pas de réseau : le chargement se replie sur `data/TSLA.csv` (déjà produit par votre `write.zoo`).
Autre actif : changer `CONFIG$ticker` et `CONFIG$fichier_cache` dans `R/00_utils.R`.

## Architecture

| Fichier | Rôle |
|:--|:--|
| `R/00_utils.R` | Configuration, packages, conventions, helpers d'affichage |
| `R/01_donnees.R` | Téléchargement + repli cache, rendements, contrôle qualité, découpage |
| `R/02_exploration.R` | Statistiques, ADF/PP/KPSS, ARCH-LM, corrélogrammes, décroissance ACF |
| `R/03_memoire_longue.R` | Périodogramme, GPH, local Whittle, R/S de Lo, Hurst (R/S, DFA, variance agrégée), d glissant |
| `R/04_farima_coeur.R` | Autocovariances exactes, Durbin-Levinson, simulation Davies-Harte, reparamétrisation |
| `R/05_farima_estimation.R` | Whittle, MLE exact (moyenne GLS profilée), écarts-types, test d = 0, résidus exacts |
| `R/06_selection.R` | Grille (p, q), critères, FARIMA contre ARMA |
| `R/07_diagnostics.R` | Portmanteau fort/faible (Lobato), ACF robuste, adéquation spectrale, stabilité, feuille de contrôle |
| `R/08_prevision.R` | Prévision AR(∞), intervalles, prix, vérification par prédicteur exact |
| `R/09_backtest.R` | Backtest en pseudo temps réel, DM (HLN), Clark-West, R² hors échantillon |
| `R/10_extensions.R` | Bootstrap par blocs, FARIMA-GARCH/FIGARCH (rugarch), VaR, Kupiec |
| `tests/test_modules.R` | Tests de correction (références indépendantes) |
| `tests/reference_python/` | Implémentation de référence en Python des formules clés (vérifiée numériquement) |

## Convention du modèle

Φ(B)(1−B)^d (X_t − μ) = Θ(B) ε_t, avec Φ(B) = 1 − Σφ_iB^i et Θ(B) = 1 + Σθ_jB^j (convention de `stats::arima`).
Ajouter un estimateur = fournir une fonction qui renvoie un objet de classe `farima_fit`.

## État de vérification

Le code R de ce projet n'a pas pu être exécuté lors de sa rédaction (R absent de l'environnement de travail).
Les formules ont été validées numériquement en Python (`tests/reference_python/`) : autocovariances contre
intégrale spectrale, vraisemblance de Durbin-Levinson contre calcul gaussien direct, moyenne GLS, simulation
Davies-Harte, écarts-types de Whittle, prévision contre prédicteur exact. La syntaxe R a été contrôlée
statiquement (équilibre des délimiteurs). **Lancer `tests/test_modules.R` en premier** : c'est la validation
du code R lui-même.
