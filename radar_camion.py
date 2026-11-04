"""Radar sur camion : que subit chaque carte électronique quand le camion roule sur une piste ?

Chemin des vibrations :  piste -> châssis -> plots (silentblocs) -> caisson -> cartes

Question du bureau d'études : QUELS PLOTS choisir pour que toutes les cartes tiennent,
sans que le caisson tape contre la structure ?
Pour chaque choix de plots, le programme donne : le Grms de chaque étage,
le débattement du caisson et la marge de sécurité de chaque carte.

Valeurs ILLUSTRATIVES, inspirées d'un radar sur camion : ce ne sont pas des données Thales.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

G = 9.81  # gravité [m/s²]

# ---------------- 1. DONNÉES (illustratives) ----------------
masse_caisson = 40.0  # [kg] caisson sans les cartes

# cartes : [nom, masse kg (avec composants), longueur mm, largeur mm, épaisseur mm, amortissement,
#           plus gros composant : longueur mm et constante C de Steinberg selon son boîtier]
# C : 0.75 axial/CMS, 1.0 DIP, 1.26 deux rangées de broches, 1.75 BGA, 2.25 LCCC céramique sans broches
cartes = [["Carte alimentation", 0.60, 160.0, 100.0, 2.0, 0.03, 30.0, 1.0],  # gros composants à broches
          ["Carte traitement", 0.35, 233.0, 160.0, 2.4, 0.03, 35.0, 1.75],  # processeur en BGA
          ["Carte RF", 0.25, 100.0, 80.0, 2.0, 0.03, 20.0, 2.25]]  # boîtier hyperfréquence céramique
E_FR4, POISSON_FR4 = 18.6e9, 0.12  # matériau des cartes (verre-époxy FR4) [Pa], [-]

# choix de plots : [nom, fréquence Hz du caisson posé sur ses plots, amortissement]
options_plots = [["Fixation rigide", 150.0, 0.05],
                 ["Plots standard 25 Hz", 25.0, 0.10],
                 ["Plots souples 10 Hz", 10.0, 0.10],
                 ["Plots souples 10 Hz amortis", 10.0, 0.25]]

# piste : l'énergie d'un camion est surtout entre 5 et 100 Hz (suspensions, roues, moteur)
profil_f = [5, 10, 40, 100, 500]  # [Hz]
profil_dsp = [0.015, 0.04, 0.04, 0.004, 0.0004]  # [g²/Hz]

jeu_mm = 10.0  # place disponible autour du caisson [mm]
f = np.geomspace(1, 1000, 20000)  # fréquences de calcul [Hz]


# ---------------- 2. CALCUL ----------------
def dsp_profil(f):
    """DSP du châssis [g²/Hz] : droites en log-log entre les points du profil, nulle en dehors"""
    dsp = 10 ** np.interp(np.log10(f), np.log10(profil_f), np.log10(profil_dsp))
    dsp[(f < profil_f[0]) | (f > profil_f[-1])] = 0
    return dsp


def rms(dsp, f):
    """valeur efficace (RMS) = racine de l'aire sous la DSP"""
    return np.sqrt(np.trapezoid(dsp, f))


def frequence_carte(carte):
    """1re fréquence propre d'une carte vue comme une plaque mince en appui simple sur ses 4 bords
    (Techniques de l'Ingénieur E3954, tableau 1, cas SSSS) : fn = π/2 √(D/ρh) (1/a² + 1/b²)"""
    nom, masse, a, b, h = carte[:5]
    a, b, h = a / 1000, b / 1000, h / 1000  # en mètres
    D = E_FR4 * h**3 / (12 * (1 - POISSON_FR4**2))  # rigidité en flexion de la plaque [N.m]
    masse_surfacique = masse / (a * b)  # ρh : masse des composants répartie sur la carte [kg/m²]
    return np.pi / 2 * np.sqrt(D / masse_surfacique) * (1 / a**2 + 1 / b**2)


def admissible_steinberg_mm(carte):
    """flexion 3σ admissible (Steinberg) : Z = 0.00022 B / (C h r √L) en pouces, pour ~20 millions d'alternances.
    r = 1 : composant supposé au centre de la carte (le cas le plus sévère)"""
    B, h, L, C = carte[2] / 25.4, carte[4] / 25.4, carte[6] / 25.4, carte[7]
    return 0.00022 * B / (C * h * 1.0 * np.sqrt(L)) * 25.4


def reponse(f_plots, a_plots, f):
    """renvoie : fréquences propres [Hz], déplacements z et transmissibilité T
    (colonne 0 = caisson, puis une colonne par carte)"""
    masses = [masse_caisson] + [c[1] for c in cartes]
    n = len(masses)
    M = np.diag(masses)
    K = np.zeros((n, n))
    C = np.zeros((n, n))

    # plots : ressort entre le châssis et le caisson (k = m ω², c = 2 ζ √(k m))
    k = sum(masses) * (2 * np.pi * f_plots) ** 2
    K[0, 0] = k
    C[0, 0] = 2 * a_plots * np.sqrt(k * sum(masses))

    # chaque carte : un ressort entre le caisson (0) et la carte (i)
    for i in range(1, n):
        m, a = cartes[i - 1][1], cartes[i - 1][5]
        k = m * (2 * np.pi * frequence_carte(cartes[i - 1])) ** 2
        c = 2 * a * np.sqrt(k * m)
        K[0, 0] += k
        K[i, i] += k
        K[0, i] -= k
        K[i, 0] -= k
        C[0, 0] += c
        C[i, i] += c
        C[0, i] -= c
        C[i, 0] -= c

    modes = np.sort(np.sqrt(np.linalg.eigvals(np.linalg.inv(M) @ K).real)) / (2 * np.pi)
    w = 2 * np.pi * f
    # à chaque fréquence : (K - ω²M + iωC) z = -M·1·a   (z = déplacement par rapport au châssis)
    z = np.array([np.linalg.solve(K - wi**2 * M + 1j * wi * C, -M @ np.ones(n)) for wi in w])
    T = 1 - w[:, None] ** 2 * z  # accélération de la pièce / accélération du châssis
    return modes, z, T


def verdict(ms):
    return "OK" if ms >= 0 else "NON CONFORME"


# ---------------- 3. PROGRAMME PRINCIPAL ----------------
if __name__ == "__main__":
    dsp = dsp_profil(f)
    grms_chassis = rms(dsp, f)
    lignes = ["RADAR SUR CAMION : CHOIX DES PLOTS DU CAISSON (valeurs illustratives)", "=" * 70, "",
              f"Entrée : châssis sur piste, {profil_f[0]}-{profil_f[-1]} Hz, {grms_chassis:.2f} Grms",
              f"Caisson : {masse_caisson} kg + {len(cartes)} cartes ; jeu disponible {jeu_mm} mm", "",
              "Cartes (plaques FR4 en appui simple, fréquence calculée ; admissible de Steinberg) :"]
    for carte in cartes:
        lignes.append(f"    {carte[0]:20s} {carte[2]:.0f}x{carte[3]:.0f}x{carte[4]} mm, {carte[1]} kg"
                      f" -> {frequence_carte(carte):4.0f} Hz ; C = {carte[7]} -> admissible {admissible_steinberg_mm(carte):.3f} mm")
    lignes.append("")
    resultats = []  # pour chaque option : [T, grms, pire marge, critère limitant]

    for nom_plots, f_plots, a_plots in options_plots:
        modes, z, T = reponse(f_plots, a_plots, f)
        grms = [rms(abs(T[:, i]) ** 2 * dsp, f) for i in range(len(modes))]
        debattement = 3 * rms(abs(z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # 3σ [mm]
        marges = {"Débattement caisson": jeu_mm / debattement - 1}  # MS = limite / calculé - 1
        lignes += [f"--- {nom_plots} ({f_plots:.0f} Hz, amortissement {a_plots:.0%})",
                   f"    Modes : {', '.join(f'{m:.0f}' for m in modes)} Hz",
                   f"    Caisson : {grms[0]:.2f} Grms, débattement 3σ {debattement:.1f} mm"
                   f"  -> MS = {marges['Débattement caisson']:+.2f}  {verdict(marges['Débattement caisson'])}"]
        for i, carte in enumerate(cartes, start=1):
            flexion = 3 * rms(abs(z[:, i] - z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # carte / caisson, 3σ [mm]
            marges[carte[0]] = admissible_steinberg_mm(carte) / flexion - 1
            lignes.append(f"    {carte[0]:28s} {grms[i]:5.2f} Grms  flexion 3σ {flexion:.3f} mm  "
                          f"MS = {marges[carte[0]]:+6.2f}  {verdict(marges[carte[0]])}")
        limitant = min(marges, key=marges.get)
        lignes += [f"    PIRE MARGE : {marges[limitant]:+.2f} ({limitant})", ""]
        resultats.append([T, grms, marges[limitant], limitant])

    lignes.append("CONCLUSION")
    for (nom_plots, _, _), (T, grms, pire, limitant) in zip(options_plots, resultats):
        lignes.append(f"  {nom_plots:28s} {verdict(pire):12s} limité par : {limitant:28s}"
                      f" pire carte = {max(grms[1:]) / grms_chassis:.1f} x châssis")
    conformes = [i for i in range(len(options_plots)) if resultats[i][2] >= 0]
    meilleure_marge = max(conformes, key=lambda i: resultats[i][2])
    mieux_protegees = min(conformes, key=lambda i: max(resultats[i][1][1:]))
    lignes += [f"  -> Plus grande marge          : {options_plots[meilleure_marge][0]}",
               f"  -> Cartes les mieux protégées : {options_plots[mieux_protegees][0]}",
               "  Le choix final est un compromis : protéger les cartes (plots souples et amortis)",
               "  sans dépasser le jeu disponible autour du caisson (débattement).",
               "  Limites : 1 axe vertical, chaque carte réduite à son 1er mode de plaque (bords en appui simple),",
               "  composants supposés au centre (r = 1), antenne et mât non modélisés, profil de piste illustratif.",
               "  À confirmer par éléments finis et mesures sur camion."]
    rapport = "\n".join(lignes)
    print(rapport)
    dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
    os.makedirs(dossier, exist_ok=True)
    with open(os.path.join(dossier, "radar_camion_rapport.txt"), "w", encoding="utf-8") as fichier:
        fichier.write(rapport)

    # figure de synthèse
    noms = [o[0] for o in options_plots]
    figure, axes = plt.subplots(2, 2, figsize=(15, 9.5))
    dans_profil = dsp > 0
    axes[0, 1].loglog(f[dans_profil], dsp[dans_profil], "k", lw=2.5, label=f"Châssis {grms_chassis:.2f} Grms")
    for nom, (T, grms, pire, limitant) in zip(noms, resultats):
        axes[0, 0].loglog(f, abs(T[:, 0]), label=nom)
        axes[0, 1].loglog(f[dans_profil], (abs(T[:, 0]) ** 2 * dsp)[dans_profil],
                          label=f"Caisson, {nom} : {grms[0]:.2f} Grms")
    axes[0, 0].axhline(1, color="grey", ls="--")
    axes[0, 0].axvspan(10, 40, color="tab:red", alpha=0.08, label="Zone d'énergie de la piste")
    axes[0, 0].set(title="Transmissibilité châssis -> caisson", xlabel="Fréquence [Hz]", xlim=(2, 1000))
    axes[0, 1].set(title="DSP d'accélération du caisson [g²/Hz]", xlabel="Fréquence [Hz]")
    for axe in axes[0]:
        axe.legend(fontsize=8)
        axe.grid(True, which="both", alpha=0.3)

    etages = ["Châssis", "Caisson"] + [c[0] for c in cartes]
    for numero, (nom, (T, grms, pire, limitant)) in enumerate(zip(noms, resultats)):
        axes[1, 0].bar(np.arange(len(etages)) + (numero - 1.5) * 0.2, [grms_chassis] + grms, width=0.2, label=nom)
    axes[1, 0].set_xticks(np.arange(len(etages)), etages, fontsize=9)
    axes[1, 0].set(title="Niveau reçu à chaque étage de la chaîne", ylabel="Grms")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, axis="y", alpha=0.3)

    pires = [r[2] for r in resultats]
    axes[1, 1].barh(noms, pires, color=["tab:green" if p >= 0 else "tab:red" for p in pires])
    axes[1, 1].axvline(0, color="k")
    axes[1, 1].set(title="Pire marge de sécurité de chaque choix (MS >= 0 : conforme)", xlabel="MS min")

    figure.suptitle("Radar sur camion : choix des plots du caisson électronique (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(os.path.join(dossier, "radar_camion.png"), dpi=150)
    plt.show()
