"""Visualisation 3D animÃ©e du radar sur camion : chÃ¢ssis, plots, caisson Ã©lectronique et cartes.

Le chÃ¢ssis est secouÃ© par un sinus Ã  la frÃ©quence choisie avec le curseur.
Les boutons Ã  gauche changent les plots du caisson (mÃªmes choix que radar_camion.py).
Chaque piÃ¨ce bouge avec l'amplitude et le retard calculÃ©s par radar_camion.py.

Ã€ essayer :
  - plots 25 Hz, curseur sur 25 Hz      : le caisson rebondit (rÃ©sonance, les plots amplifient)
  - plots 10 Hz amortis, curseur 25 Hz  : le caisson reste presque immobile (isolation)
  - fixation rigide, curseur sur 120 Hz : la carte alimentation (orange) flÃ©chit fortement

Les dÃ©placements sont fortement AGRANDIS et ralentis pour Ãªtre visibles.
L'antenne, la cabine et les roues sont dÃ©coratives (non modÃ©lisÃ©es).

Dans Spyder : F5. L'animation s'ouvre dans une fenÃªtre sÃ©parÃ©e (il faut PyQt6 dans .venv).
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import RadioButtons, Slider

from radar_camion import cartes, options_plots, reponse
from visualisation_3d import dessiner_boite, dessiner_ressort

try:  # dans Spyder (IPython) : fenÃªtre sÃ©parÃ©e, sinon l'animation reste figÃ©e dans le panneau Graphiques
    get_ipython().run_line_magic("matplotlib", "qt")  # noqa: F821
except NameError:
    pass

AMPLITUDE_CHASSIS = 0.12  # amplitude affichÃ©e du mouvement du chÃ¢ssis (agrandie)
AGRANDISSEMENT_CARTES = 3.0  # la flexion des cartes est en plus agrandie x3
NB_IMAGES_PAR_PERIODE = 40
COULEURS_CARTES = ["tab:orange", "tab:green", "tab:purple"]

# GÃ©omÃ©trie du dessin (unitÃ©s arbitraires, x = longueur du camion)
LONGUEUR_CHASSIS, LARGEUR_CHASSIS = 9.0, 2.4
CAISSON_X, CAISSON_LONGUEUR, CAISSON_HAUTEUR = 4.2, 3.6, 2.2
HAUTEUR_PLOTS = 1.0


def calculer_transmissibilites(indice_plots, frequence_Hz):
    """transmissibilitÃ© complexe chÃ¢ssis -> [caisson, carte 1, carte 2, carte 3] Ã  une frÃ©quence"""
    nom, f_plots, a_plots = options_plots[indice_plots]
    return reponse(f_plots, a_plots, np.array([frequence_Hz]))[2][0]


def dessiner_roue(axe, x, y, z_centre):
    angle = np.linspace(0, 2 * np.pi, 40)
    axe.plot(x + 0.55 * np.cos(angle), np.full_like(angle, y), z_centre + 0.55 * np.sin(angle), color="k", lw=3,
             zorder=0)  # dessinÃ©e sous le chÃ¢ssis


def dessiner_scene(axe, phase, indice_plots, frequence_Hz):
    """dessine le camion radar Ã  un instant donnÃ© (phase = Ï‰ t)"""
    transmissibilites = calculer_transmissibilites(indice_plots, frequence_Hz)
    oscillation = np.exp(1j * phase)
    # Ã©chelle du dessin : si une piÃ¨ce bouge beaucoup plus que le chÃ¢ssis, on rÃ©duit tout
    # pour que le dessin reste lisible (les vrais rapports sont Ã©crits dans le titre)
    echelle = AMPLITUDE_CHASSIS / max(1.0, np.abs(transmissibilites).max() / 3)
    deplacement_chassis = echelle * np.cos(phase)
    deplacements = echelle * (transmissibilites * oscillation).real  # caisson, cartes...

    axe.cla()
    # ordre d'affichage imposÃ© : les traits (ressorts) passent devant les volumes, sinon ils sont cachÃ©s
    axe.computed_zorder = False
    # roues (fixes) et chÃ¢ssis (secouÃ© par la piste)
    for x in (1.3, 6.0, 7.6):
        for y in (0.0, LARGEUR_CHASSIS):
            dessiner_roue(axe, x, y, 0.55)
    haut_chassis = 1.3 + deplacement_chassis
    dessiner_boite(axe, 0, LONGUEUR_CHASSIS, 0, LARGEUR_CHASSIS, haut_chassis - 0.3, haut_chassis, "dimgrey", 0.9)
    # cabine (dÃ©corative, liÃ©e au chÃ¢ssis)
    dessiner_boite(axe, 0.2, 2.6, 0.1, LARGEUR_CHASSIS - 0.1, haut_chassis, haut_chassis + 2.0, "lightgrey", 0.5)
    # antenne radar sur un mÃ¢t liÃ© au chÃ¢ssis (dÃ©corative)
    axe.plot([8.3, 8.3], [1.2, 1.2], [haut_chassis, haut_chassis + 3.2], color="k", lw=3)
    dessiner_boite(axe, 8.1, 8.4, -0.4, LARGEUR_CHASSIS + 0.4, haut_chassis + 3.2, haut_chassis + 5.0, "silver", 0.8)

    # caisson Ã©lectronique sur 4 plots
    fond_caisson = haut_chassis + HAUTEUR_PLOTS + (deplacements[0] - deplacement_chassis)
    for x in (CAISSON_X + 0.3, CAISSON_X + CAISSON_LONGUEUR - 0.3):
        for y in (0.4, LARGEUR_CHASSIS - 0.4):
            dessiner_ressort(axe, x, y, haut_chassis, fond_caisson)
    dessiner_boite(axe, CAISSON_X, CAISSON_X + CAISSON_LONGUEUR, 0.2, LARGEUR_CHASSIS - 0.2,
                   fond_caisson, fond_caisson + CAISSON_HAUTEUR, "tab:blue", 0.12)

    # cartes empilÃ©es dans le caisson : bords tenus, centre qui flÃ©chit (forme de cuvette)
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
    niveaux = ", ".join(f"{carte[0].replace('Carte ', '')} x{abs(t):.2f}"
                        for carte, t in zip(cartes, transmissibilites[1:]))
    axe.set_title(f"{options_plots[indice_plots][0]} - excitation {frequence_Hz:.0f} Hz\n"
                  f"amplitude / chÃ¢ssis : caisson x{abs(transmissibilites[0]):.2f} ; cartes : {niveaux}\n"
                  f"(dÃ©placements agrandis et ralentis, flexion des cartes x{AGRANDISSEMENT_CARTES:.0f} en plus)",
                  fontsize=9)


def main():
    figure = plt.figure(figsize=(12, 8))
    axe = figure.add_axes([0.22, 0.12, 0.78, 0.85], projection="3d")
    curseur = Slider(figure.add_axes([0.3, 0.04, 0.55, 0.03]), "FrÃ©quence [Hz]", 3, 300, valinit=25.0)
    noms_plots = [option[0] for option in options_plots]
    boutons = RadioButtons(figure.add_axes([0.01, 0.4, 0.2, 0.2]), noms_plots, active=1)
    legende = "\n".join(f"{carte[0]} ({carte[2]:.0f} Hz)" for carte in cartes)
    figure.text(0.01, 0.3, "Cartes, de bas en haut :\n" + legende, fontsize=9)

    def mettre_a_jour(numero_image):
        phase = 2 * np.pi * numero_image / NB_IMAGES_PAR_PERIODE
        indice_plots = noms_plots.index(boutons.value_selected)
        dessiner_scene(axe, phase, indice_plots, curseur.val)

    # garder une rÃ©fÃ©rence Ã  l'animation, sinon Python l'efface et elle s'arrÃªte
    figure.animation = FuncAnimation(figure, mettre_a_jour, frames=NB_IMAGES_PAR_PERIODE, interval=60)
    plt.show()


def enregistrer_gif(indice_plots, frequence_Hz, chemin):
    """enregistre une pÃ©riode d'animation en GIF (pour le README)"""
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
