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
| `recalage.py` | Recale le modèle sur un essai au pot vibrant par descente de gradient (essai simulé). |
| `visualisation_3d.py` | Animation 3D du caisson, des silentblocs et de la carte, avec un curseur de fréquence. |
| `maquette.cpp` | Le même calcul que `maquette.py` en C++, pour montrer les bases du C++. |
| `resultats/` | Les sorties, recréées à chaque lancement. |

Le code suit le style du cours [ceriam.net](https://ceriam.net/) :
- des noms explicites, sans variable d'une seule lettre, avec un verbe pour chaque fonction ;
- des constantes en MAJUSCULES ;
- des `@dataclass` (`Etage`, `Resultat`) pour regrouper les données ;
- des annotations de type ;
- une exception `DonneeInvalide`, levée si une donnée n'a pas de sens physique ;
- une `seed` pour que le hasard soit le même à chaque lancement.

## Lancer

**Dans Spyder** : ouvre `maquette.py`, puis appuie sur **F5**. L'interpréteur doit être `.venv\Scripts\python.exe` (*Outils → Préférences → Interpréteur Python*).

**En ligne de commande** :

```bash
.venv\Scripts\python.exe maquette.py
```

```bash
.venv\Scripts\python.exe recalage.py
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

## Visualisation 3D (`visualisation_3d.py`)

Dans Spyder, ouvre `visualisation_3d.py`, puis appuie sur **F5**. Une fenêtre séparée s'ouvre : la base vibre, les 4 silentblocs (ressorts orange) se déforment, le caisson monte et descend, et la carte (en vert) fléchit en cuvette. Le **curseur** en bas de la fenêtre change la fréquence d'excitation, de 5 à 300 Hz. Tu peux aussi faire tourner la vue à la souris.

| Fréquence | Ce qu'on voit |
|---|---|
| ~25 Hz | Résonance de la suspension : le caisson bouge 5 fois plus que la base |
| ~60 Hz | Isolation : le caisson bouge beaucoup moins que la base |
| ~200 Hz | Le caisson est presque immobile, mais la carte vibre (sa propre résonance) |

Les amplitudes et les retards viennent du calcul de `maquette.py`. Les déplacements sont agrandis et ralentis pour être visibles, et la flexion de la carte est agrandie trois fois de plus.

L'animation a besoin de PyQt6 dans `.venv` (déjà installé). Si la fenêtre reste figée, va dans *Outils → Préférences → Console IPython → Graphiques → Sortie graphique* et choisis **Automatique**.

Aperçus enregistrés : `resultats/animation_25Hz.gif`, `animation_60Hz.gif` et `animation_200Hz.gif`.

![Résonance de la suspension à 25 Hz](resultats/animation_25Hz.gif)

## Recalage modèle-essai (`recalage.py`)

La fiche technique annonce des silentblocs à 25 Hz et 10 % d'amortissement, mais le vrai produit est différent. On le mesure au pot vibrant, puis on ajuste le modèle pour qu'il colle à la mesure. C'est la mission « corrélation des modèles » de l'offre.

La méthode reprend le cours (*Forward Neural Network → Loss / Training*) :

| Notion du cours | Dans `recalage.py` |
|---|---|
| Paramètres θ | fréquence et amortissement des silentblocs (normalisés autour de 1) |
| Perte J (MSE avec le facteur ½) | écart entre la transmissibilité du modèle et celle mesurée (en log) |
| Descente de gradient θ ← θ − α∇J | `recaler()`, avec un pas α = 0,5 |
| Critère d'arrêt ‖α∇J‖ < ε | ε = 10⁻⁷ |
| Jeu d'entraînement / jeu de test | 70 % des points servent à recaler, 30 % à juger le résultat |
| Graine (`seed`) | le bruit de mesure est le même à chaque lancement |

Le gradient est calculé par **différences finies**. Avec seulement 2 paramètres, c'est suffisant : l'autodifférentiation vue dans le cours ne devient utile qu'avec beaucoup de paramètres.

**Résultat** (l'essai est simulé, avec un vrai caisson caché et 5 % de bruit) :

| | Fiche technique | Recalé | Vrai (caché) |
|---|---|---|---|
| Fréquence | 25,0 Hz | 22,2 Hz | 22,0 Hz |
| Amortissement | 10 % | 14,1 % | 14 % |

- La perte sur le jeu de test passe de 9,6·10⁻³ à 1,4·10⁻⁴.
- La pire marge de sécurité passe de +1,93 à +2,47.
- Figure : `resultats/recalage.png`.

**Suite possible, vue dans le cours** : un **PINN** (réseau de neurones informé par la physique) pourrait faire ce recalage quand le modèle physique est mal connu (élastomère non linéaire, effet de la température). Pour un système linéaire comme celui-ci, la méthode ci-dessus reste plus simple et plus précise.

## À dire en entretien

- **Isolation** : au-delà de √2 × 25 ≈ 35 Hz, les silentblocs atténuent la vibration (6 Grms deviennent 1,5 Grms). En dessous, ils l'amplifient (Q ≈ 5 à 25 Hz).
- **Amortissement, un compromis** : plus d'amortissement réduit le pic à 25 Hz, mais au-delà de 20 % environ la carte reçoit de nouveau plus d'énergie.
- **Température** : un élastomère se raidit au froid. La sensibilité montre si le produit reste conforme sur toute la plage de température.
- **Limites** : un seul axe, des masses ponctuelles, un comportement linéaire. La suite logique est un modèle éléments finis, puis un essai au pot vibrant pour recaler les fréquences et les amortissements.
