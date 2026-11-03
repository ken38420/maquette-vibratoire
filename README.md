# Maquette vibratoire

Un caisson électronique posé sur des **silentblocs**, avec une **carte électronique** à l'intérieur. Le porteur (camion, avion…) secoue le tout avec une **vibration aléatoire**. Le programme calcule ce que reçoivent le caisson et la carte, puis vérifie que tout tient.

> Toutes les valeurs sont **illustratives** : ce ne sont pas des données Thales.

```
   carte (m2)      0,3 kg   vibre seule à 200 Hz
     ↕  ressort de la carte
   caisson (m1)    25 kg    rebondit sur ses silentblocs à 25 Hz
     ↕  silentblocs
   base            secouée par le profil d'essai (6 Grms)
```

## Les fichiers

| Fichier | Rôle |
|---|---|
| `maquette.py` | Tout le modèle Python. Les données sont en haut du fichier. |
| `maquette.cpp` | Le même calcul en C++, pour montrer les bases du C++. |
| `resultats/` | Les sorties, recréées à chaque lancement. |

## Lancer

**Dans Spyder** : ouvre `maquette.py`, puis appuie sur **F5**. L'interpréteur doit être `.venv\Scripts\python.exe` (*Outils → Préférences → Interpréteur Python*).

**En ligne de commande** :

```bash
.venv\Scripts\python.exe maquette.py
```

**Version C++** :

```bash
g++ -O2 -o maquette.exe maquette.cpp
```

```bash
maquette.exe
```

Pour essayer d'autres valeurs (masse, fréquence, amortissement, profil), modifie la partie `1. DONNÉES D'ENTRÉE` en haut de `maquette.py`, puis relance.

## Ce que le programme sort (ce qu'un bureau d'études veut voir)

1. **`rapport.txt`** : une note de calcul avec six parties :
   - les hypothèses ;
   - les modes propres (fréquence, masse effective, déformée) ;
   - la réponse aléatoire (Grms, amplification Q) ;
   - les critères avec leurs marges ;
   - la sensibilité ;
   - la vérification du code.
2. **`synthese.png`** : quatre graphiques :
   - la transmissibilité ;
   - les DSP d'entrée et de sortie ;
   - les marges du cas nominal ;
   - la pire marge de chaque cas de sensibilité.
3. **`courbes.csv`** : les courbes en colonnes, à ouvrir dans Excel (séparateur `;`).

**Les trois critères vérifiés.** Pour chacun, on calcule une marge de sécurité **MS = admissible / calculé − 1**. Le critère est conforme si MS ≥ 0.

| Critère | Question posée |
|---|---|
| Séparation des modes | La carte vibre-t-elle au moins 2 fois plus haut que la suspension (règle de l'octave) ? |
| Débattement des silentblocs (3σ) | Le caisson risque-t-il de taper contre la structure (jeu de 5 mm) ? |
| Déplacement de la carte (3σ, Steinberg) | Les soudures des composants vont-elles casser par fatigue ? |

**Sensibilité.** Le programme refait le calcul dans quatre autres cas :
- silentblocs plus mous (au chaud, −30 % de raideur) ;
- silentblocs plus raides (au froid, +50 % de raideur) ;
- amortissement de 5 % ;
- amortissement de 20 %.

Il garde ensuite le **pire cas**.

**Vérification du code.** Avant de faire confiance aux résultats, le programme vérifie un cas dont on connaît la réponse exacte : un oscillateur simple sous bruit blanc (équation de Miles). L'écart obtenu est de 3·10⁻⁵.

## Résultats (cas nominal)

| | Valeur |
|---|---|
| Modes | 25,0 Hz (le caisson sur ses silentblocs) et 201,2 Hz (la carte) |
| Niveaux | entrée 6,06 Grms → caisson 1,47 Grms → carte 1,56 Grms |
| Débattement 3σ | 1,71 mm pour un jeu de 5 mm (MS = +1,93) |
| Carte 3σ | 0,029 mm pour 0,563 mm admissibles (MS = +18,4) |
| Pire cas | amortissement faible (5 %) : MS = +0,98, toujours conforme |

## À dire en entretien

- **Isolation** : au-delà de √2 × 25 ≈ 35 Hz, les silentblocs atténuent la vibration (6 Grms deviennent 1,5 Grms). En dessous, ils l'amplifient (Q ≈ 5 à 25 Hz).
- **Amortissement, un compromis** : plus d'amortissement réduit le pic à 25 Hz, mais au-delà de 20 % environ la carte reçoit de nouveau plus d'énergie.
- **Température** : un élastomère se raidit au froid. La sensibilité montre si le produit reste conforme sur toute la plage de température.
- **Limites** : un seul axe, des masses ponctuelles, un comportement linéaire. La suite logique est un modèle éléments finis, puis un essai au pot vibrant pour recaler les fréquences et les amortissements.
