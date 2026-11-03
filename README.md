# Maquette vibratoire — caisson sur silentblocs + carte électronique

Modèle à 2 degrés de liberté pour illustrer, en 5 minutes d'entretien, la démarche du « modèle vibratoire 3D d'un produit en caisson » de l'offre Thales : modes propres, transmissibilité, réponse à une DSP d'entrée, effet de l'amortissement, critère de tenue d'une carte.

> **Valeurs illustratives.** Masses, fréquences, amortissements et profil de DSP sont des hypothèses de démo, pas des données Thales ni un profil AECTP 401 réel. Tout est dans `config.json`.

## Lancer

Le projet est géré avec [uv](https://docs.astral.sh/uv/), comme dans le cours ceriam.net. Depuis le dossier `maquette_vibratoire` :

```bash
uv sync
```

```bash
uv run maquette-vibratoire config.json
```

`uv sync` crée l'environnement `.venv` avec les versions exactes de `uv.lock`. Sorties : résumé dans la console, `figures/synthese.png` et `figures/resultats.json` (paramètres, résultats et versions des bibliothèques).

> Sur ce PC, le réseau intercepte le HTTPS : si `uv sync` échoue sur un certificat, lancer d'abord `$env:UV_SYSTEM_CERTS='1'` dans PowerShell.

## Tests

```bash
uv run pytest -v
```

```bash
uv run ruff check src tests
```

24 tests pytest, chacun comparant le code à un résultat connu par ailleurs :

| Test | Référence |
|---|---|
| Paramètres invalides (masse ≤ 0, ζ hors ]0, 1[, profil incohérent, CLI sans argument) | Le code doit s'arrêter avec un message clair |
| Modes propres | Formule fermée du 2 ddl |
| 1 ddl sous bruit blanc | Solution exacte `Grms² = S₀ π fₙ (1+4ζ²)/(4ζ)` |
| 1 ddl à la résonance | `|T| = √(1+4ζ²)/(2ζ)`, pour 4 valeurs de ζ |
| 1 ddl à √2·fₙ | `|T| = 1` quel que soit ζ (début de l'isolation) |
| 2 ddl à f → 0 | `T = 1` (tout suit la base) |
| DSP trapèze | +3 dB/octave : une octave sous le palier, la DSP vaut niveau × 10^(-0,3) |
| RMS, Grms d'entrée | Intégrales calculées analytiquement |
| Steinberg | B = h = L = 1 pouce donne 0,00022 pouce |
| Convergence | Grille divisée par deux, écart < 1e-3 |
| **Non-régression** | L'étude complète redonne, à 1e-9 près, les résultats de la version validée (`tests/donnees/reference.json`) |

Le test de non-régression tourne dans un dossier temporaire (fixture `tmp_path`) et n'écrase donc pas les figures du projet.

## Le modèle

```
   x2 ── carte (m2)        k2, c2 : liaison carte / caisson  (f ≈ 200 Hz, ζ = 3 %)
   x1 ── caisson (m1)      k1, c1 : silentblocs              (f ≈ 25 Hz,  ζ = 10 %)
   y  ── base (porteur)    accélération imposée = DSP d'entrée
```

- Équation en coordonnées relatives : `M z'' + C z' + K z = -M·1·ÿ`
- FRF calculée directement : `Z(ω) = (K - ω²M + iωC)⁻¹ (-M·1)`, transmissibilité `T = 1 - ω² Z`
- Réponse aléatoire : `S_sortie = |T|² · S_entrée`, puis `Grms = √∫S df`
- Tenue de la carte : déplacement relatif carte/caisson à 3σ comparé au critère de Steinberg

## Résultats (config actuelle)

| Grandeur | Valeur |
|---|---|
| Mode 1 (caisson sur silentblocs) | 25,0 Hz, caisson et carte en phase |
| Mode 2 (carte) | 201,2 Hz, carte en opposition |
| Entrée | 6,06 Grms |
| Caisson / carte | 1,47 / 1,56 Grms |
| Déplacement carte 3σ / admissible Steinberg | 0,029 / 0,563 mm (marge ≈ 19) |

## Validations intégrées (le script s'arrête si l'une échoue)

1. **Modes** : `eigh` comparé à la formule fermée du 2 ddl (écart < 1e-9).
2. **1 ddl sous bruit blanc** : le même code est comparé à la solution exacte `Grms² = S₀ π fₙ (1+4ζ²)/(4ζ)`, c'est-à-dire l'équation de Miles corrigée (écart ≈ 4e-6).
3. **Convergence** : le Grms de la carte ne change pas quand on divise la grille fréquentielle par deux (écart ≈ 5e-5).
4. **Limite basse fréquence** : T(f→0) = 1, la masse suit la base.

## Ce qu'il faut dire en entretien

- **Isolation** : au-delà de √2 × 25 ≈ 35 Hz, les silentblocs atténuent. L'entrée passe de 6 Grms à 1,5 Grms sur le caisson. En dessous, ils amplifient (×5 à 25 Hz). Il faut donc placer la fréquence de suspension sous le début du profil.
- **Amortissement, compromis** : plus d'amortissement réduit le pic à 25 Hz, mais au-delà de ζ ≈ 0,2 la carte reçoit *plus* d'énergie, car l'isolation haute fréquence se dégrade (T ≈ 2ζ fₙ/f). Ce n'est pas « plus il y en a, mieux c'est ».
- **Découplage** : la carte (200 Hz) est bien séparée de la suspension (25 Hz), les deux modes interagissent peu. Règle courante : un rapport d'octave (×2) au minimum entre deux modes successifs.
- **Lien avec le vrai modèle 3D** : la démarche est identique avec un modèle éléments finis (Ansys, Abaqus, Nastran). On passe de la CAO à la géométrie simplifiée, puis au maillage, aux conditions aux limites, à l'analyse modale, à la réponse aléatoire, et enfin au **recalage** sur un essai au pot vibrant (fréquences, amortissement).
- **Lien qualification** : le profil d'entrée viendrait de l'AECTP 400 (méthode 401) adapté au porteur. La maquette permet de chercher quels essais mutualiser et quelles marges vérifier avant l'essai.

## Limites (à assumer)

- 1 direction, 2 masses ponctuelles : pas de modes de flexion du caisson, pas de couplage entre axes.
- Amortissement visqueux constant. Un élastomère réel dépend de la fréquence, de l'amplitude et de la température.
- Le critère de Steinberg est empirique (composants standards). Pour un vrai dossier, on fait un calcul de fatigue (Miner, Dirlik) sur les contraintes du modèle EF.

## Version C++ (`cpp/`)

Même modèle que la version Python, écrit en C++17 sans aucune bibliothèque externe, pour montrer la maîtrise du C++ demandée par l'offre.

| Partie | Comment c'est fait en C++ |
|---|---|
| Modes propres | Formule fermée de `det(K - λM) = 0`. Les racines sont calculées sous une forme stable qui évite l'annulation numérique. |
| FRF | `std::complex<double>`, système 2×2 résolu par la règle de Cramer |
| DSP + intégration | Profil trapèze en dB/octave, grille logarithmique, méthode des trapèzes |
| Validation | 1 ddl sous bruit blanc comparé à la solution exacte (écart 3,5e-6) |
| Erreurs | Exceptions (`std::invalid_argument`), message `ERREUR : …` et code de sortie 1 |

Compiler et lancer (GCC 16.2, WinLibs), depuis le dossier `cpp` :

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -static -o modele_vibratoire.exe modele_vibratoire.cpp
```

```bash
modele_vibratoire.exe 25 0.3 25 0.10 200 0.03
```

Les arguments sont, dans l'ordre : `m1_kg m2_kg f_susp_Hz zeta_susp f_carte_Hz zeta_carte`. Le profil de DSP et la grille sont écrits en dur, avec les mêmes valeurs que `config.json`.

**Comparaison avec Python** (même config) : modes 25,00 et 201,22 Hz, transmissibilités max 5,1 et 5,2, accélérations 6,06 / 1,47 / 1,56 Grms, déplacement 3σ 0,029 mm. Les résultats sont identiques à la précision affichée.

## Fichiers

```
maquette_vibratoire/
├── pyproject.toml            dépendances figées (numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2), pytest, ruff
├── uv.lock                   versions exactes de toutes les dépendances
├── config.json               tous les paramètres (aucune valeur par défaut cachée dans le code)
├── src/maquette_vibratoire/
│   └── modele.py             modèle, validations, figures (annoté en types)
├── tests/
│   ├── conftest.py           fixtures : config, référence, dossier temporaire
│   ├── test_modele.py        24 tests
│   └── donnees/reference.json
├── figures/                  sorties
└── cpp/modele_vibratoire.cpp version C++ (calcul seul, sans figure)
```
