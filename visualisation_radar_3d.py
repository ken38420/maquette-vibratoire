"""Visualisation 3D animée du radar sur camion : châssis, plots, caisson électronique et cartes.

Le châssis est secoué par un sinus à la fréquence choisie avec le curseur.
Les boutons à gauche changent les plots du caisson (mêmes choix que radar_camion.py).
Chaque pièce bouge avec l'amplitude et le retard calculés par radar_camion.py.

À essayer :
  - plots 25 Hz, curseur sur 25 Hz      : le caisson rebondit (résonance, les plots amplifient)
  - plots 10 Hz amortis, curseur 25 Hz  : le caisson reste presque immobile (isolation)
  - fixation rigide, curseur sur 120 Hz : la carte alimentation (orange) fléchit fortement

Les déplacements sont fortement AGRANDIS et ralentis pour être visibles.
L'antenne, la cabine et les roues sont décoratives (non modélisées).

Dans Spyder : F5. L'animation s'ouvre dans une fenêtre séparée (il faut PyQt6 dans .venv).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import RadioButtons, Slider

from radar_camion import CARTES, OPTIONS_PLOTS, construire_matrices
from visualisation_3d import dessiner_boite, dessiner_ressort

try:  # dans Spyder (IPython) : fenêtre séparée, sinon l'animation reste figée dans le panneau Graphiques
    get_ipython().run_line_magic("matplotlib", "qt")  # noqa: F821
except NameError:
    pass

AMPLITUDE_CHASSIS = 0.12  # amplitude affichée du mouvement du châssis (agrandie)
AGRANDISSEMENT_CARTES = 3.0  # la flexion des cartes est en plus agrandie x3
NB_IMAGES_PAR_PERIODE = 40
COULEURS_CARTES = ["tab:orange", "tab:green", "tab:purple"]

# Géométrie du dessin (unités arbitraires, x = longueur du camion)
LONGUEUR_CHASSIS, LARGEUR_CHASSIS = 9.0, 2.4
CAISSON_X, CAISSON_LONGUEUR, CAISSON_HAUTEUR = 4.2, 3.6, 2.2
HAUTEUR_PLOTS = 1.0


def calculer_transmissibilites(indice_plots: int, frequence_Hz: float) -> np.ndarray:
    """transmissibilité complexe châssis -> [caisson, carte 1, carte 2, carte 3] à une fréquence"""
    matrice_masse, matrice_raideur, matrice_amort = construire_matrices(OPTIONS_PLOTS[indice_plots])
    omega = 2 * np.pi * frequence_Hz
    deplacements = np.linalg.solve(matrice_raideur - omega**2 * matrice_masse + 1j * omega * matrice_amort,
                                   -matrice_masse @ np.ones(len(matrice_masse)))
    return 1 - omega**2 * deplacements


def dessiner_roue(axe, x, y, z_centre) -> None:
    angle = np.linspace(0, 2 * np.pi, 40)
    axe.plot(x + 0.55 * np.cos(angle), np.full_like(angle, y), z_centre + 0.55 * np.sin(angle), color="k", lw=3,
             zorder=0)  # dessinée sous le châssis


def dessiner_scene(axe, phase: float, indice_plots: int, frequence_Hz: float) -> None:
    """dessine le camion radar à un instant donné (phase = ω t)"""
    transmissibilites = calculer_transmissibilites(indice_plots, frequence_Hz)
    oscillation = np.exp(1j * phase)
    # échelle du dessin : si une pièce bouge beaucoup plus que le châssis, on réduit tout
    # pour que le dessin reste lisible (les vrais rapports sont écrits dans le titre)
    echelle = AMPLITUDE_CHASSIS / max(1.0, np.abs(transmissibilites).max() / 3)
    deplacement_chassis = echelle * np.cos(phase)
    deplacements = echelle * (transmissibilites * oscillation).real  # caisson, cartes...

    axe.cla()
    # ordre d'affichage imposé : les traits (ressorts) passent devant les volumes, sinon ils sont cachés
    axe.computed_zorder = False
    # roues (fixes) et châssis (secoué par la piste)
    for x in (1.3, 6.0, 7.6):
        for y in (0.0, LARGEUR_CHASSIS):
            dessiner_roue(axe, x, y, 0.55)
    haut_chassis = 1.3 + deplacement_chassis
    dessiner_boite(axe, 0, LONGUEUR_CHASSIS, 0, LARGEUR_CHASSIS, haut_chassis - 0.3, haut_chassis, "dimgrey", 0.9)
    # cabine (décorative, liée au châssis)
    dessiner_boite(axe, 0.2, 2.6, 0.1, LARGEUR_CHASSIS - 0.1, haut_chassis, haut_chassis + 2.0, "lightgrey", 0.5)
    # antenne radar sur un mât lié au châssis (décorative)
    axe.plot([8.3, 8.3], [1.2, 1.2], [haut_chassis, haut_chassis + 3.2], color="k", lw=3)
    dessiner_boite(axe, 8.1, 8.4, -0.4, LARGEUR_CHASSIS + 0.4, haut_chassis + 3.2, haut_chassis + 5.0, "silver", 0.8)

    # caisson électronique sur 4 plots
    fond_caisson = haut_chassis + HAUTEUR_PLOTS + (deplacements[0] - deplacement_chassis)
    for x in (CAISSON_X + 0.3, CAISSON_X + CAISSON_LONGUEUR - 0.3):
        for y in (0.4, LARGEUR_CHASSIS - 0.4):
            dessiner_ressort(axe, x, y, haut_chassis, fond_caisson)
    dessiner_boite(axe, CAISSON_X, CAISSON_X + CAISSON_LONGUEUR, 0.2, LARGEUR_CHASSIS - 0.2,
                   fond_caisson, fond_caisson + CAISSON_HAUTEUR, "tab:blue", 0.12)

    # cartes empilées dans le caisson : bords tenus, centre qui fléchit (forme de cuvette)
    u, v = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 10))
    forme = np.sin(np.pi * u) * np.sin(np.pi * v)
    for numero, couleur in enumerate(COULEURS_CARTES):
        hauteur_bords = fond_caisson + 0.5 + 0.55 * numero
        flexion = AGRANDISSEMENT_CARTES * (deplacements[numero + 1] - deplacements[0])
        axe.plot_surface(CAISSON_X + 0.4 + u * (CAISSON_LONGUEUR - 0.8), 0.5 + v * (LARGEUR_CHASSIS - 1.0),
                         hauteur_bords + flexion * forme, color=couleur, alpha=0.9)

    axe.set(xlim=(-0.5, LONGUEUR_CHASSIS + 0.5), ylim=(-1.5, LARGEUR_CHASSIS + 1.5), zlim=(0, 7))
    axe.set_box_aspect((LONGUEUR_CHASSIS + 1, LARGEUR_CHASSIS + 3, 7))
    axe.view_init(elev=15, azim=-70)
    axe.set_axis_off()
    niveaux = ", ".join(f"{carte.nom.replace('Carte ', '')} x{abs(t):.2f}"
                        for carte, t in zip(CARTES, transmissibilites[1:], strict=True))
    axe.set_title(f"{OPTIONS_PLOTS[indice_plots].nom} - excitation {frequence_Hz:.0f} Hz\n"
                  f"amplitude / châssis : caisson x{abs(transmissibilites[0]):.2f} ; cartes : {niveaux}\n"
                  f"(déplacements agrandis et ralentis, flexion des cartes x{AGRANDISSEMENT_CARTES:.0f} en plus)",
                  fontsize=9)


def main() -> None:
    figure = plt.figure(figsize=(12, 8))
    axe = figure.add_axes([0.22, 0.12, 0.78, 0.85], projection="3d")
    curseur = Slider(figure.add_axes([0.3, 0.04, 0.55, 0.03]), "Fréquence [Hz]", 3, 300, valinit=25.0)
    boutons = RadioButtons(figure.add_axes([0.01, 0.4, 0.2, 0.2]), [plots.nom for plots in OPTIONS_PLOTS], active=1)
    legende = "\n".join(f"{carte.nom} ({carte.frequence_Hz:.0f} Hz)" for carte in CARTES)
    figure.text(0.01, 0.3, "Cartes, de bas en haut :\n" + legende, fontsize=9)

    def mettre_a_jour(numero_image: int) -> None:
        phase = 2 * np.pi * numero_image / NB_IMAGES_PAR_PERIODE
        indice_plots = [plots.nom for plots in OPTIONS_PLOTS].index(boutons.value_selected)
        dessiner_scene(axe, phase, indice_plots, curseur.val)

    # garder une référence à l'animation, sinon Python l'efface et elle s'arrête
    figure.animation = FuncAnimation(figure, mettre_a_jour, frames=NB_IMAGES_PAR_PERIODE, interval=60)
    plt.show()


def enregistrer_gif(indice_plots: int, frequence_Hz: float, chemin: str) -> None:
    """enregistre une période d'animation en GIF (pour le README)"""
    figure = plt.figure(figsize=(9, 6.5))
    axe = figure.add_subplot(projection="3d")
    animation = FuncAnimation(
        figure,
        lambda numero: dessiner_scene(axe, 2 * np.pi * numero / NB_IMAGES_PAR_PERIODE, indice_plots, frequence_Hz),
        frames=NB_IMAGES_PAR_PERIODE,
    )
    animation.save(chemin, writer="pillow", fps=20)
    plt.close(figure)


if __name__ == "__main__":
    main()
