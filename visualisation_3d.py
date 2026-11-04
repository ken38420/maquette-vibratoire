"""Visualisation 3D animée de la maquette : base, silentblocs, caisson et carte.

La base est secouée par un sinus à la fréquence choisie avec le curseur.
Chaque pièce bouge avec l'amplitude et le retard calculés par maquette.py
(transmissibilité complexe T) :
  - vers 25 Hz  : le caisson rebondit sur ses silentblocs (résonance, il bouge plus que la base) ;
  - vers 60 Hz  : les silentblocs isolent (le caisson bouge beaucoup moins que la base) ;
  - vers 200 Hz : le caisson est presque immobile, mais la carte vibre (résonance de la carte).

Les déplacements sont fortement AGRANDIS et ralentis pour être visibles.

Dans Spyder : F5. L'animation s'ouvre dans une fenêtre séparée (il faut PyQt6 dans .venv).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from maquette import caisson, carte, reponse

try:  # dans Spyder (IPython) : fenêtre séparée, sinon l'animation reste figée dans le panneau Graphiques
    get_ipython().run_line_magic("matplotlib", "qt")  # noqa: F821
except NameError:
    pass

AMPLITUDE_BASE = 0.25  # amplitude affichée du mouvement de la base (unité du dessin, agrandie)
AGRANDISSEMENT_CARTE = 3.0  # la flexion de la carte est en plus agrandie x3 pour être visible
NB_IMAGES_PAR_PERIODE = 40  # ralenti : une période d'oscillation = 40 images
FREQUENCE_INITIALE_HZ = 25.0

# Géométrie du dessin (unités arbitraires)
LARGEUR, PROFONDEUR = 4.0, 3.0  # caisson
HAUTEUR_CAISSON = 2.0
HAUTEUR_SILENTBLOCS = 1.0
HAUTEUR_CARTE = 0.8  # hauteur de la carte au-dessus du fond du caisson


def calculer_transmissibilites(frequence_Hz):
    """transmissibilité complexe (amplitude et retard) du caisson et de la carte à une fréquence"""
    T = reponse(caisson, carte, np.array([frequence_Hz]))[2][0]
    return T[0], T[1]


def dessiner_boite(axe, x0, x1, y0, y1, z0, z1, couleur, transparence):
    """parallélépipède défini par ses coins"""
    coins = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                      [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    faces = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]]
    axe.add_collection3d(Poly3DCollection([coins[face] for face in faces], facecolor=couleur,
                                          edgecolor="k", linewidth=0.5, alpha=transparence))


def dessiner_ressort(axe, x, y, z_bas, z_haut):
    """silentbloc dessiné comme un ressort hélicoïdal entre z_bas et z_haut"""
    angle = np.linspace(0, 5 * 2 * np.pi, 200)  # 5 spires : assez espacées pour voir le ressort s'étirer
    rayon = 0.25
    hauteurs = np.linspace(z_bas, z_haut, angle.size)
    axe.plot(x + rayon * np.cos(angle), y + rayon * np.sin(angle), hauteurs, color="black", linewidth=5)  # contour
    axe.plot(x + rayon * np.cos(angle), y + rayon * np.sin(angle), hauteurs, color="red", linewidth=3)


def dessiner_scene(axe, phase, frequence_Hz):
    """dessine toutes les pièces à un instant donné (phase = ω t)"""
    transmissibilite_caisson, transmissibilite_carte = calculer_transmissibilites(frequence_Hz)
    oscillation = np.exp(1j * phase)
    # déplacement vertical de chaque pièce = partie réelle de (amplitude complexe × e^{iωt})
    deplacement_base = AMPLITUDE_BASE * np.cos(phase)
    deplacement_caisson = AMPLITUDE_BASE * (transmissibilite_caisson * oscillation).real
    deplacement_carte = AMPLITUDE_BASE * (transmissibilite_carte * oscillation).real

    axe.cla()
    # ordre d'affichage imposé : les traits (ressorts) passent devant les volumes, sinon ils sont cachés
    axe.computed_zorder = False
    # base (le porteur)
    dessiner_boite(axe, -0.5, LARGEUR + 0.5, -0.5, PROFONDEUR + 0.5,
                   deplacement_base - 0.2, deplacement_base, "grey", 0.35)
    # 4 silentblocs aux coins
    fond_caisson = HAUTEUR_SILENTBLOCS + deplacement_caisson
    for x, y in [(0.3, 0.3), (LARGEUR - 0.3, 0.3), (0.3, PROFONDEUR - 0.3), (LARGEUR - 0.3, PROFONDEUR - 0.3)]:
        dessiner_ressort(axe, x, y, deplacement_base, fond_caisson)
    # caisson transparent
    dessiner_boite(axe, 0, LARGEUR, 0, PROFONDEUR, fond_caisson, fond_caisson + HAUTEUR_CAISSON, "tab:blue", 0.12)

    # carte : bords fixés au caisson, centre qui bouge avec la masse "carte" -> forme de cuvette
    u, v = np.meshgrid(np.linspace(0, 1, 12), np.linspace(0, 1, 12))
    forme = np.sin(np.pi * u) * np.sin(np.pi * v)  # 0 sur les bords, 1 au centre
    hauteur_bords = fond_caisson + HAUTEUR_CARTE
    flexion = AGRANDISSEMENT_CARTE * (deplacement_carte - deplacement_caisson)  # centre par rapport aux bords
    axe.plot_surface(0.5 + u * (LARGEUR - 1), 0.5 + v * (PROFONDEUR - 1), hauteur_bords + flexion * forme,
                     color="tab:green", alpha=0.9)

    axe.set(xlim=(-0.5, LARGEUR + 0.5), ylim=(-0.5, PROFONDEUR + 0.5), zlim=(-0.6, 4.2))
    axe.set_box_aspect((LARGEUR + 1, PROFONDEUR + 1, 4.8))
    axe.view_init(elev=12, azim=-60)  # vue presque de côté pour bien voir le mouvement vertical
    axe.set_axis_off()
    axe.set_title(f"Excitation {frequence_Hz:.0f} Hz   |   amplitude / base : caisson x{abs(transmissibilite_caisson):.2f},"
                  f" carte x{abs(transmissibilite_carte):.2f}\n(déplacements agrandis et ralentis,"
                  f" flexion de la carte x{AGRANDISSEMENT_CARTE:.0f} en plus)")


def main():
    figure = plt.figure(figsize=(9, 8))
    axe = figure.add_subplot(projection="3d")
    figure.subplots_adjust(bottom=0.12)
    curseur = Slider(figure.add_axes([0.2, 0.04, 0.6, 0.03]), "Fréquence [Hz]", 5, 300,
                     valinit=FREQUENCE_INITIALE_HZ)

    def mettre_a_jour(numero_image):
        phase = 2 * np.pi * numero_image / NB_IMAGES_PAR_PERIODE
        dessiner_scene(axe, phase, curseur.val)

    # garder une référence à l'animation, sinon Python l'efface et elle s'arrête
    figure.animation = FuncAnimation(figure, mettre_a_jour, frames=NB_IMAGES_PAR_PERIODE, interval=50)
    plt.show()


def enregistrer_gif(frequence_Hz, chemin):
    """enregistre une période d'animation en GIF (pour le README)"""
    figure = plt.figure(figsize=(7, 6.5))
    axe = figure.add_subplot(projection="3d")
    animation = FuncAnimation(figure, lambda numero: dessiner_scene(axe, 2 * np.pi * numero / NB_IMAGES_PAR_PERIODE,
                                                                     frequence_Hz),
                              frames=NB_IMAGES_PAR_PERIODE)
    animation.save(chemin, writer="pillow", fps=20)
    plt.close(figure)


if __name__ == "__main__":
    main()
