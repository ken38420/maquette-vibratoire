# Fiche d'entretien : stage qualification environnementale, Thales

## À savoir avant d'entrer

- **Le site.** Thales **Fleury-les-Aubrais**, près d'Orléans : le principal site de support et de services pour les radars et les systèmes d'armes. Ne pas confondre avec **Ymare**, près de Rouen, qui produit les radars Ground Master.
- **Les radars.** Plus de 300 Ground Master vendus dans plus de 40 pays. Chiffres trouvés en ligne, à revérifier la veille.
- **Ce que veut le responsable.** L'offre parle de « rationalisation », de « standardiser le référentiel » et de « réduire les délais de qualification ». Il attend donc surtout de la **méthode de qualification**, pas seulement de la simulation.
- **Question à lui poser.** L'offre parle d'un « nouvel outil de maintenance » : de quoi s'agit-il ?

## Ordre de la démo (5 minutes)

1. `radar_camion.png` : le problème, les 4 choix de plots et le compromis.
2. `visualisation_radar_3d.py` : plots 25 Hz à 25 Hz, puis fixation rigide à 135 Hz.
3. `matrice_essais.png` : passer de 20 essais à 5 grâce à l'enveloppe et aux preuves existantes.
4. `recalage.png` : comment on valide le modèle avec un essai.

Garder les PNG et les GIF ouverts à l'avance, au cas où Spyder poserait problème.

## Les normes, en une phrase chacune

| Norme | Contenu |
|---|---|
| AECTP-100 | Stratégie : le **profil de vie** du matériel (transport, stockage, utilisation), d'où l'on tire les essais |
| AECTP-300 | Climatique : chaud, froid, chaleur humide (306), brouillard salin (309), sable, pluie… |
| AECTP-400 | Mécanique : vibration (401), chocs… |
| AECTP-500 | CEM : émissions et susceptibilité, renvoie à MIL-STD-461 |
| Feu | **Pas dans l'AECTP-300** : il faut une exigence spécifique (réaction au feu des matériaux, fumées) |

Les AECTP ne fixent pas de niveau universel : on **personnalise** l'essai à partir du profil de vie. L'équivalent américain est la MIL-STD-810 (la méthode 514 pour les vibrations). Pour un véhicule à roues, on trouve de l'ordre de 1 à 2 Grms et environ 40 min par axe : ordre de grandeur à vérifier dans la norme.

Les numéros de méthode marqués ici sont ceux de mes recherches. Les vérifier avant de les citer avec assurance.

## Questions probables

1. **« Comment construiriez-vous un référentiel de qualification commun ? »**
   Je recense les exigences de chaque projet et je les regroupe par environnement. J'en prends l'enveloppe, je définis des niveaux standard, puis je relie chaque exigence à une preuve : un essai, une analyse ou une similitude avec un produit déjà qualifié. C'est ce que fait `matrice_essais.py`.

2. **« Comment passe-t-on d'un profil de vie à un essai ? »**
   On part du profil de vie (AECTP-100) : les scénarios, les durées, les distances et les spectres mesurés. On en tire l'enveloppe, puis on **compresse la durée** par équivalence de dommage. Avec la loi de Basquin, T_essai = T_vie · (G_vie / G_essai)^b, avec b entre 4 et 8 environ. Il ne faut pas trop augmenter le niveau, sinon on crée des défaillances qui n'existent pas en vrai.

3. **« Comment mutualiser les essais ? »**
   On fait un essai enveloppe pour plusieurs projets. On réutilise les preuves d'un produit déjà qualifié quand la conception est la même. On enchaîne les essais sur une seule éprouvette : la mécanique d'abord, le climatique ensuite, le brouillard salin en dernier, et le feu à part car il est destructif.

4. **« Quels sont les critères de réussite d'un essai ? »**
   Le produit fonctionne avant, pendant et après l'essai, sans dégradation visible. Ces critères sont écrits **avant** l'essai, dans le plan d'essais, et reliés à l'exigence. On vérifie aussi que le pot vibrant suit bien la consigne (tolérance de ±3 dB sur la DSP).

5. **« Comment validez-vous un modèle ? »**
   D'abord, je vérifie le code sur un cas dont je connais la réponse exacte (l'oscillateur sous bruit blanc). Ensuite, un essai de balayage sinus à bas niveau. Puis je recale les fréquences et les amortissements. Enfin, je contrôle l'écart sur un **autre** essai que celui qui a servi au recalage.

6. **« Pourquoi des plots ne sont-ils pas toujours bénéfiques ? »**
   En dessous de √2 × leur fréquence, ils amplifient. Sur un camion, l'énergie est entre 10 et 40 Hz, donc des plots à 25 Hz amplifient. Il faut viser une fréquence bien plus basse, et vérifier que le caisson ne tape pas.

7. **« Qu'est-ce que la règle de l'octave ? »**
   Chaque mode doit être au moins 2 fois plus haut que le mode qui le porte. Sinon les deux se couplent. Dans ma maquette, sans plots, la carte traitement reçoit 6,7 fois le niveau du châssis.

8. **« Quelles sont les limites de votre maquette ? »**
   Un seul axe, des masses ponctuelles, un comportement linéaire et des données inventées. La suite serait un modèle éléments finis du caisson (j'ai utilisé Abaqus à Polimi), puis un recalage sur un essai au pot vibrant.

## Questions pièges sur la physique

- **Pourquoi la carte reçoit-elle plus de Grms que le caisson ?**
  Presque tout vient de la zone sous 100 Hz : la carte suit le mouvement à 25 Hz du caisson, comme un bloc rigide. Sa propre résonance à 200 Hz n'ajoute que peu.

- **Pourquoi |T| = 1 à √2 × fn ?**
  Pour un système à un degré de liberté, la transmissibilité vaut 1 à r = √2, quel que soit l'amortissement. Au-delà, on isole.

- **Pourquoi multiplier la flexion par 1,62 ?**
  Une masse sur un ressort donne le déplacement « moyen » de la carte. Au centre d'une plaque en appui sur ses 4 bords, la flexion est 16/π² ≈ 1,62 fois plus grande.

- **Pourquoi 3σ pour Steinberg et 4σ pour le débattement ?**
  Steinberg a établi son critère en 3σ, il faut donc garder cette convention. Pour savoir si le caisson tape, il faut le plus grand pic de l'essai : avec une loi de Rayleigh, sur une heure d'essai à 25 Hz (environ 90 000 cycles), il atteint environ 4 à 4,5σ.

- **Miles ou solution exacte ?**
  Miles, Grms ≈ √(π/2·fn·Q·S), est une approximation valable pour un amortissement faible et une DSP plate. Le code utilise la solution exacte avec le facteur (1+4ζ²). Pour le caisson, Miles donne 1,57 Grms contre 1,47 calculé, car la DSP n'est pas plate autour de 25 Hz.

- **Amortissement visqueux ou hystérétique pour un élastomère ?**
  Plutôt hystérétique : k(1+iη), avec η ≈ 2ζ. Au-delà de la résonance, le modèle visqueux décroît moins vite. Il est donc prudent (conservatif) pour l'isolation.

- **Que change la température ?**
  Au froid, l'élastomère se raidit et s'amortit souvent davantage. Raideur et amortissement varient donc ensemble, alors que le code ne fait varier que la raideur. Ce serait une amélioration à faire.

- **Pourquoi l'amortissement est-il un compromis ?**
  Plus d'amortissement réduit le pic à la résonance, mais transmet plus au-delà. Dans la maquette, le Grms de la carte est minimal vers 20 % d'amortissement, puis il remonte.

## C++ (l'offre le demande)

`maquette.cpp` refait le calcul de `maquette.py` et donne les mêmes chiffres. À savoir expliquer :
- `std::complex<double>` pour les nombres complexes ;
- la règle de Cramer pour résoudre le système 2 × 2 ;
- la méthode des trapèzes pour intégrer la DSP ;
- la compilation avec `g++ -O2 -o maquette.exe maquette.cpp`.
