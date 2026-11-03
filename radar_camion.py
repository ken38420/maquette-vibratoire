"""Radar monté sur camion : que subit chaque carte électronique quand le camion roule sur une piste ?

Chaîne de transmission des vibrations :

    piste -> châssis du camion -> plots (silentblocs) -> caisson électronique -> cartes

  - le châssis est secoué par la piste : c'est l'ENTRÉE (une DSP d'accélération, en g²/Hz) ;
  - le caisson est posé sur des plots, qui filtrent (ou amplifient !) les vibrations ;
  - chaque carte est fixée dans le caisson et a sa propre fréquence de résonance.

Question du bureau d'études : QUELS PLOTS choisir pour que toutes les cartes tiennent,
sans que le caisson tape contre la structure ?

Le programme compare plusieurs choix de plots et sort, pour chacun :
Grms de chaque étage, débattement du caisson, marge de sécurité de chaque carte.

Valeurs ILLUSTRATIVES, inspirées d'un radar sur camion (type Ground Master) :
ce ne sont pas des données Thales.
"""

import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

type Vecteur = npt.NDArray[np.float64]

GRAVITE = 9.81  # [m/s²]
MM_PAR_POUCE = 25.4


@dataclass
class Carte:
    nom: str
    masse_kg: float
    frequence_Hz: float  # fréquence propre de la carte fixée sur un caisson rigide
    amortissement: float
    bord_mm: float  # longueur du bord de la carte (pour Steinberg)
    epaisseur_mm: float
    composant_mm: float  # longueur du plus gros composant


@dataclass
class Plots:
    nom: str
    frequence_Hz: float  # fréquence du caisson complet posé sur ses plots
    amortissement: float


# ============================================================
# 1. DONNÉES D'ENTRÉE (illustratives)
# ============================================================
MASSE_CAISSON_KG = 40.0  # caisson + châssis interne, sans les cartes

CARTES = [
    Carte("Carte alimentation", 0.60, 120.0, 0.03, 160.0, 2.0, 30.0),
    Carte("Carte traitement", 0.35, 250.0, 0.03, 233.0, 1.6, 35.0),  # traitement du signal
    Carte("Carte RF", 0.25, 400.0, 0.03, 100.0, 1.6, 20.0),  # hyperfréquence
]

OPTIONS_PLOTS = [
    Plots("Fixation rigide", 150.0, 0.05),
    Plots("Plots standard 25 Hz", 25.0, 0.10),
    Plots("Plots souples 10 Hz", 10.0, 0.10),
    Plots("Plots souples 10 Hz amortis", 10.0, 0.25),
]

# Profil de piste : l'énergie d'un camion est surtout entre 5 et 100 Hz
# (suspensions, roues, moteur). Points (fréquence [Hz], DSP [g²/Hz]), droites en log-log.
PROFIL_FREQUENCES_HZ = [5, 10, 40, 100, 500]
PROFIL_DSP_G2_HZ = [0.015, 0.04, 0.04, 0.004, 0.0004]

JEU_MM = 10.0  # débattement disponible autour du caisson [mm]
FREQUENCES_HZ = np.geomspace(1, 1000, 20000)


# ============================================================
# 2. CALCUL
# ============================================================
def calculer_dsp_entree(frequences: Vecteur) -> Vecteur:
    """DSP du châssis [g²/Hz] : interpolation log-log du profil, nulle en dehors"""
    log_dsp = np.interp(np.log10(frequences), np.log10(PROFIL_FREQUENCES_HZ), np.log10(PROFIL_DSP_G2_HZ))
    dsp = 10**log_dsp
    dsp[(frequences < PROFIL_FREQUENCES_HZ[0]) | (frequences > PROFIL_FREQUENCES_HZ[-1])] = 0.0
    return dsp


def calculer_rms(dsp: Vecteur, frequences: Vecteur) -> float:
    """valeur efficace (RMS) = racine de l'aire sous la DSP"""
    return float(np.sqrt(np.trapezoid(dsp, frequences)))


def construire_matrices(plots: Plots) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """matrices M, K, C : degré de liberté 0 = caisson, 1, 2, 3... = cartes"""
    nb_ddl = 1 + len(CARTES)
    masse_suspendue = MASSE_CAISSON_KG + sum(carte.masse_kg for carte in CARTES)
    matrice_masse = np.diag([MASSE_CAISSON_KG] + [carte.masse_kg for carte in CARTES])
    matrice_raideur = np.zeros((nb_ddl, nb_ddl))
    matrice_amort = np.zeros((nb_ddl, nb_ddl))

    # plots : entre le châssis (imposé) et le caisson -> k = m ω², c = 2 ζ √(k m)
    raideur_plots = masse_suspendue * (2 * np.pi * plots.frequence_Hz) ** 2
    matrice_raideur[0, 0] += raideur_plots
    matrice_amort[0, 0] += 2 * plots.amortissement * np.sqrt(raideur_plots * masse_suspendue)

    # chaque carte : un ressort entre le caisson (0) et la carte (indice)
    for indice, carte in enumerate(CARTES, start=1):
        raideur = carte.masse_kg * (2 * np.pi * carte.frequence_Hz) ** 2
        amort = 2 * carte.amortissement * np.sqrt(raideur * carte.masse_kg)
        for matrice, valeur in ((matrice_raideur, raideur), (matrice_amort, amort)):
            matrice[0, 0] += valeur
            matrice[indice, indice] += valeur
            matrice[0, indice] -= valeur
            matrice[indice, 0] -= valeur
    return matrice_masse, matrice_raideur, matrice_amort


def calculer_reponse(plots: Plots, frequences: Vecteur, dsp_entree: Vecteur) -> dict:
    """réponse de toute la chaîne pour un choix de plots"""
    matrice_masse, matrice_raideur, matrice_amort = construire_matrices(plots)
    un = np.ones(len(matrice_masse))

    # modes propres
    omega_carre = np.linalg.eigvals(np.linalg.inv(matrice_masse) @ matrice_raideur).real
    frequences_propres = np.sort(np.sqrt(omega_carre)) / (2 * np.pi)

    # réponse en fréquence : (K - ω²M + iωC) z = -M·1·a_châssis
    deplacements = np.array(
        [
            np.linalg.solve(matrice_raideur - omega**2 * matrice_masse + 1j * omega * matrice_amort, -matrice_masse @ un)
            for omega in 2 * np.pi * frequences
        ]
    )
    transmissibilite = 1 - (2 * np.pi * frequences[:, None]) ** 2 * deplacements
    dsp_sortie = np.abs(transmissibilite) ** 2 * dsp_entree[:, None]  # colonne 0 caisson, puis cartes

    dsp_debattement = np.abs(deplacements[:, 0]) ** 2 * dsp_entree * GRAVITE**2
    flexions_3sigma_mm = [
        3 * calculer_rms(np.abs(deplacements[:, indice] - deplacements[:, 0]) ** 2 * dsp_entree * GRAVITE**2,
                         frequences) * 1000
        for indice in range(1, len(un))
    ]
    return {
        "frequences_propres_Hz": frequences_propres,
        "transmissibilite_caisson": np.abs(transmissibilite[:, 0]),
        "dsp_sortie": dsp_sortie,
        "grms": [calculer_rms(dsp_sortie[:, indice], frequences) for indice in range(len(un))],
        "debattement_3sigma_mm": 3 * calculer_rms(dsp_debattement, frequences) * 1000,
        "flexions_3sigma_mm": flexions_3sigma_mm,
    }


def calculer_admissible_steinberg_mm(carte: Carte) -> float:
    """critère de Steinberg (formule en pouces) : au-dessous, ~10 millions de cycles avant rupture des soudures"""
    bord, epaisseur, composant = (x / MM_PAR_POUCE for x in (carte.bord_mm, carte.epaisseur_mm, carte.composant_mm))
    return 0.00022 * bord / (epaisseur * np.sqrt(composant)) * MM_PAR_POUCE


def calculer_marges(reponse: dict) -> dict[str, float]:
    """marges de sécurité MS = admissible / calculé - 1 (conforme si MS >= 0)"""
    marges = {"Débattement caisson": JEU_MM / reponse["debattement_3sigma_mm"] - 1}
    for carte, flexion in zip(CARTES, reponse["flexions_3sigma_mm"], strict=True):
        marges[carte.nom] = calculer_admissible_steinberg_mm(carte) / flexion - 1
    return marges


# ============================================================
# 3. SORTIES
# ============================================================
def conforme(marge: float) -> str:
    return "OK" if marge >= 0 else "NON CONFORME"


def ecrire_rapport(grms_chassis: float, reponses: list[dict], marges: list[dict]) -> str:
    lignes = [
        "RADAR SUR CAMION : CHOIX DES PLOTS DU CAISSON (valeurs illustratives)", "=" * 70, "",
        f"Entrée : châssis du camion sur piste, {PROFIL_FREQUENCES_HZ[0]}-{PROFIL_FREQUENCES_HZ[-1]} Hz, "
        f"{grms_chassis:.2f} Grms",
        f"Caisson : {MASSE_CAISSON_KG} kg + {len(CARTES)} cartes ; jeu disponible {JEU_MM} mm", "",
    ]
    for plots, reponse, marge in zip(OPTIONS_PLOTS, reponses, marges, strict=True):
        lignes += [
            f"--- {plots.nom} ({plots.frequence_Hz:.0f} Hz, amortissement {plots.amortissement:.0%})",
            f"    Modes : {', '.join(f'{f:.0f}' for f in reponse['frequences_propres_Hz'])} Hz",
            f"    Caisson : {reponse['grms'][0]:.2f} Grms, débattement 3σ {reponse['debattement_3sigma_mm']:.1f} mm"
            f"  -> MS = {marge['Débattement caisson']:+.2f}  {conforme(marge['Débattement caisson'])}",
        ]
        for indice, carte in enumerate(CARTES, start=1):
            lignes.append(f"    {carte.nom:28s} {reponse['grms'][indice]:5.2f} Grms  "
                          f"flexion 3σ {reponse['flexions_3sigma_mm'][indice - 1]:.3f} mm  "
                          f"MS = {marge[carte.nom]:+6.2f}  {conforme(marge[carte.nom])}")
        critere_limitant = min(marge, key=marge.get)
        lignes.append(f"    PIRE MARGE : {marge[critere_limitant]:+.2f} ({critere_limitant})")
        lignes.append("")

    lignes.append("CONCLUSION")
    for plots, reponse, marge in zip(OPTIONS_PLOTS, reponses, marges, strict=True):
        critere_limitant = min(marge, key=marge.get)
        rapport_cartes = max(reponse["grms"][1:]) / grms_chassis  # pire carte / châssis
        lignes.append(f"  {plots.nom:28s} {conforme(marge[critere_limitant]):12s} limité par : {critere_limitant:28s}"
                      f" pire carte = {rapport_cartes:.1f} x châssis")
    conformes = [indice for indice in range(len(OPTIONS_PLOTS)) if min(marges[indice].values()) >= 0]
    if conformes:
        plus_grande_marge = max(conformes, key=lambda indice: min(marges[indice].values()))
        cartes_protegees = min(conformes, key=lambda indice: max(reponses[indice]["grms"][1:]))
        lignes += [
            f"  -> Plus grande marge             : {OPTIONS_PLOTS[plus_grande_marge].nom}",
            f"  -> Cartes les mieux protégées    : {OPTIONS_PLOTS[cartes_protegees].nom}",
            "  Le choix final est un compromis : protéger les cartes (plots souples et amortis)",
            "  sans dépasser le jeu disponible autour du caisson (débattement).",
        ]
    lignes += [
        "  Limites : 1 axe vertical, cartes = masses ponctuelles, antenne et mât non modélisés,",
        "  profil de piste illustratif. À confirmer par un modèle éléments finis et des mesures sur camion.",
    ]
    return "\n".join(lignes)


def tracer_synthese(dsp_entree, grms_chassis, reponses, marges, chemin: str) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(15, 9.5))
    frequences = FREQUENCES_HZ

    axe = axes[0, 0]
    for plots, reponse in zip(OPTIONS_PLOTS, reponses, strict=True):
        axe.loglog(frequences, reponse["transmissibilite_caisson"], label=plots.nom)
    axe.axhline(1, color="grey", ls="--")
    axe.axvspan(PROFIL_FREQUENCES_HZ[1], PROFIL_FREQUENCES_HZ[2], color="tab:red", alpha=0.08,
                label="Zone d'énergie de la piste")
    axe.set(title="Transmissibilité châssis -> caisson", xlabel="Fréquence [Hz]", xlim=(2, 1000))
    axe.legend(fontsize=8)
    axe.grid(True, which="both", alpha=0.3)

    axe = axes[0, 1]
    dans_profil = dsp_entree > 0
    axe.loglog(frequences[dans_profil], dsp_entree[dans_profil], "k", lw=2.5, label=f"Châssis {grms_chassis:.2f} Grms")
    for plots, reponse in zip(OPTIONS_PLOTS, reponses, strict=True):
        axe.loglog(frequences[dans_profil], reponse["dsp_sortie"][dans_profil, 0],
                   label=f"Caisson, {plots.nom} : {reponse['grms'][0]:.2f} Grms")
    axe.set(title="DSP d'accélération du caisson [g²/Hz]", xlabel="Fréquence [Hz]")
    axe.legend(fontsize=8)
    axe.grid(True, which="both", alpha=0.3)

    axe = axes[1, 0]
    etages = ["Châssis", "Caisson"] + [carte.nom for carte in CARTES]
    largeur = 0.8 / len(OPTIONS_PLOTS)
    for numero, (plots, reponse) in enumerate(zip(OPTIONS_PLOTS, reponses, strict=True)):
        positions = np.arange(len(etages)) + (numero - (len(OPTIONS_PLOTS) - 1) / 2) * largeur
        axe.bar(positions, [grms_chassis] + reponse["grms"], width=largeur, label=plots.nom)
    axe.set_xticks(np.arange(len(etages)), etages, fontsize=9)
    axe.set(title="Niveau reçu à chaque étage de la chaîne", ylabel="Grms")
    axe.legend(fontsize=8)
    axe.grid(True, axis="y", alpha=0.3)

    axe = axes[1, 1]
    pires_marges = [min(marge.values()) for marge in marges]
    axe.barh([plots.nom for plots in OPTIONS_PLOTS], pires_marges,
             color=["tab:green" if marge >= 0 else "tab:red" for marge in pires_marges])
    axe.axvline(0, color="k")
    axe.set(title="Pire marge de sécurité de chaque choix (MS >= 0 : conforme)", xlabel="MS min")

    figure.suptitle("Radar sur camion : choix des plots du caisson électronique (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(chemin, dpi=150)


def main() -> None:
    dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
    os.makedirs(dossier, exist_ok=True)

    dsp_entree = calculer_dsp_entree(FREQUENCES_HZ)
    grms_chassis = calculer_rms(dsp_entree, FREQUENCES_HZ)
    reponses = [calculer_reponse(plots, FREQUENCES_HZ, dsp_entree) for plots in OPTIONS_PLOTS]
    marges = [calculer_marges(reponse) for reponse in reponses]

    rapport = ecrire_rapport(grms_chassis, reponses, marges)
    print(rapport)
    with open(os.path.join(dossier, "radar_camion_rapport.txt"), "w", encoding="utf-8") as fichier:
        fichier.write(rapport)
    tracer_synthese(dsp_entree, grms_chassis, reponses, marges, os.path.join(dossier, "radar_camion.png"))
    plt.show()


if __name__ == "__main__":
    main()
