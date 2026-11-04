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
| `radar_camion.py` | **Exemple pour l'entretien** : radar sur camion, chaîne châssis → plots → caisson → 3 cartes, comparaison de 4 choix de plots. |
| `visualisation_radar_3d.py` | Animation 3D du camion radar, avec un curseur de fréquence et des boutons pour choisir les plots. |
| `maquette.cpp` | Le même calcul que `maquette.py` en C++, pour montrer les bases du C++. |
| `resultats/` | Les sorties, recréées à chaque lancement. |

Le code est volontairement **le plus simple possible** :
- seulement `numpy` et `matplotlib` ;
- les données sont de simples listes en haut du fichier, par exemple `caisson = [masse, fréquence, amortissement]` ;
- une seule fonction de calcul, `reponse()`, qui construit les matrices M, K, C et résout le système à chaque fréquence ;
- le reste du fichier se lit de haut en bas : calcul, note de calcul, figure.

Le calcul tient en une ligne : `(K − ω²M + iωC) z = −M·1·a`, résolue à chaque fréquence.

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

Pour essayer d'autres valeurs (masse, fréquence, amortissement, profil), modifie la partie `1. DONNÉES` en haut de `maquette.py`, puis relance.

## Ce que le programme sort (ce qu'un bureau d'études veut voir)

1. **`rapport.txt`** : une note de calcul avec six parties :
   - les hypothèses ;
   - les fréquences propres ;
   - les niveaux reçus (Grms, amplification maximale) ;
   - les critères avec leurs marges ;
   - la sensibilité ;
   - la vérification du code.
2. **`synthese.png`** : quatre graphiques :
   - la transmissibilité ;
   - les DSP d'entrée et de sortie ;
   - les marges du cas nominal ;
   - la pire marge de chaque cas de sensibilité.

**Les trois critères vérifiés.** Pour chacun, on calcule une marge de sécurité **MS = admissible / calculé − 1** (pour la séparation des modes, qui est un minimum : calculé / minimum − 1). Le critère est conforme si MS ≥ 0.

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

## Exemple pour l'entretien : radar sur camion (`radar_camion.py`)

Prenons un radar monté sur camion, du type Ground Master. L'antenne est sur un mât. L'électronique est rangée dans des caissons posés sur des **plots** (silentblocs). Quand le camion roule sur une piste, les vibrations suivent une chaîne :

```
piste → châssis du camion → plots → caisson électronique → cartes (alimentation, traitement, RF)
        (ENTRÉE : profil de piste)  (filtrent ou amplifient)   (chacune a sa propre résonance)
```

**La question du bureau d'études :** quels plots choisir pour que toutes les cartes tiennent, sans que le caisson tape contre la structure ?

> Valeurs **illustratives** : profil de piste, masses, fréquences et dimensions des cartes sont inventés pour la démonstration. Ce ne sont pas des données du Ground Master.

**Ce que montre le modèle** (profil de piste à 1,6 Grms, énergie surtout entre 10 et 40 Hz) :

| Choix de plots | Caisson | Pire carte | Débattement 3σ (jeu 10 mm) | Verdict |
|---|---|---|---|---|
| Fixation rigide | 2,7 Grms | **8,4 Grms** (alimentation) | 0,1 mm | **Non conforme** : le mode du caisson (≈ 150 Hz) se couple avec la carte alimentation (120 Hz) |
| Plots standard 25 Hz | **2,8 Grms** | 2,9 Grms | 3,3 mm | Conforme, mais les plots **amplifient** : 25 Hz tombe en plein dans l'énergie de la piste |
| Plots souples 10 Hz | 1,7 Grms | 1,7 Grms | **12 mm** | **Non conforme** : le caisson tape contre la structure |
| Plots souples 10 Hz amortis (25 %) | 1,1 Grms | 1,1 Grms | 6,9 mm | Conforme, et ce sont les cartes les mieux protégées |

**Les leçons à retenir :**
1. **Des plots ne sont pas « bons » ou « mauvais » dans l'absolu.** Des plots à 25 Hz isolent bien sur un avion, dont l'énergie est à haute fréquence, mais ils amplifient sur un camion. Il faut toujours partir du **profil de vie** (AECTP-100).
2. **Deux modes proches se couplent.** Sans plots, le mode du caisson tombe sur la carte alimentation, qui reçoit 5 fois le niveau du châssis. C'est la règle de l'octave, vue en pratique.
3. **Le choix final est un compromis.** Des plots souples protègent les cartes mais demandent du jeu autour du caisson. L'amortissement permet de réduire le débattement.

Sorties : `resultats/radar_camion_rapport.txt` (note de calcul) et `resultats/radar_camion.png` (4 graphiques : transmissibilité, DSP du caisson, Grms à chaque étage, pire marge).

**Animation 3D** : `visualisation_radar_3d.py` (F5 dans Spyder). Le curseur change la fréquence et les boutons à gauche changent les plots. À essayer :

| Réglage | Ce qu'on voit |
|---|---|
| Plots 25 Hz, curseur sur 25 Hz | le caisson rebondit sur ses plots (×5) |
| Plots 10 Hz amortis, curseur sur 25 Hz | le caisson reste presque immobile (×0,3) |
| Fixation rigide, curseur sur 120 Hz | la carte alimentation (orange) fléchit très fortement |

Aperçus : `resultats/radar_plots25Hz_25Hz.gif`, `radar_plots10Hz_amortis_25Hz.gif` et `radar_rigide_120Hz.gif`.

![Radar sur camion, plots 25 Hz excités à 25 Hz](resultats/radar_plots25Hz_25Hz.gif)

## Visualisation 3D (`visualisation_3d.py`)

Dans Spyder, ouvre `visualisation_3d.py`, puis appuie sur **F5**. Une fenêtre séparée s'ouvre : la base vibre, les 4 silentblocs (ressorts rouges) se déforment, le caisson monte et descend, et la carte (en vert) fléchit en cuvette. Le **curseur** en bas de la fenêtre change la fréquence d'excitation, de 5 à 300 Hz. Tu peux aussi faire tourner la vue à la souris.

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
| Descente de gradient θ ← θ − α∇J | la boucle `for` du fichier, avec un pas α = 0,5 |
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
