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
# C : 0.75 composants axiaux, 1.0 DIP (à broches), 1.26 deux rangées de broches, 1.75 BGA, 2.25 LCCC céramique
cartes = [["Carte alimentation", 0.60, 160.0, 100.0, 2.0, 0.03, 30.0, 1.0],  # gros composants à broches
          ["Carte traitement", 0.35, 233.0, 160.0, 2.4, 0.03, 35.0, 1.75],  # processeur en BGA
          ["Carte RF", 0.25, 100.0, 80.0, 2.0, 0.03, 20.0, 2.25]]  # boîtier hyperfréquence céramique
E_FR4, POISSON_FR4 = 18.6e9, 0.12  # matériau des cartes (verre-époxy FR4) [Pa], [-]

# choix de plots : [nom, fréquence Hz du caisson posé sur ses plots, amortissement]
options_plots = [["Fixation rigide", 150.0, 0.05],
                 ["Plots standard 25 Hz", 25.0, 0.10],
                 ["Plots souples 10 Hz", 10.0, 0.10],
                 ["Plots souples 10 Hz amortis", 10.0, 0.25]]

# piste : le plateau (le plus d'énergie) est entre 10 et 40 Hz (suspensions, roues), puis décroît jusqu'à 500 Hz
profil_f = [5, 10, 40, 100, 500]  # [Hz]
profil_dsp = [0.015, 0.04, 0.04, 0.004, 0.0004]  # [g²/Hz]

jeu_mm = 10.0  # place disponible autour du caisson [mm]

# chaque carte est réduite à UNE masse sur UN ressort : cela donne son déplacement "moyen" (modal).
# Au centre d'une plaque en appui sur ses 4 bords, la flexion est 16/π² ≈ 1,62 fois plus grande.
FACTEUR_CENTRE = 16 / np.pi**2
f = np.geomspace(1, 1000, 20000)  # fréquences de calcul [Hz]


# ---------------- 2. CALCUL ----------------
# dsp_profil, rms, reponse et verdict : mêmes idées que dans maquette.py, avec 3 cartes au lieu d'une
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
    nom, masse, longueur, largeur, h = carte[:5]
    a, b, h = longueur / 1000, largeur / 1000, h / 1000  # côtés a, b et épaisseur h, en mètres
    D = E_FR4 * h**3 / (12 * (1 - POISSON_FR4**2))  # rigidité en flexion de la plaque [N.m]
    masse_surfacique = masse / (a * b)  # ρh : masse des composants répartie sur la carte [kg/m²]
    return np.pi / 2 * np.sqrt(D / masse_surfacique) * (1 / a**2 + 1 / b**2)


def admissible_steinberg_mm(carte):
    """flexion 3σ admissible (Steinberg) : Z = 0.00022 B / (C h r √L) en pouces, pour ~20 millions d'alternances.
    r = 1 : composant supposé au centre de la carte (le cas le plus sévère)"""
    B, h, L, C = carte[2] / 25.4, carte[4] / 25.4, carte[6] / 25.4, carte[7]
    return 0.00022 * B / (C * h * 1.0 * np.sqrt(L)) * 25.4


def reponse(f_plots, amort_plots, f):
    """renvoie : fréquences propres [Hz], déplacements z et transmissibilité T
    (colonne 0 = caisson, puis une colonne par carte)"""
    masses = [masse_caisson] + [c[1] for c in cartes]
    n = len(masses)
    M = np.diag(masses)
    K = np.zeros((n, n))
    C = np.zeros((n, n))

    # plots : ressort entre le châssis et le caisson (k = m ω², c = 2 ζ √(k m))
    k_plots = sum(masses) * (2 * np.pi * f_plots) ** 2
    K[0, 0] = k_plots
    C[0, 0] = 2 * amort_plots * np.sqrt(k_plots * sum(masses))

    # chaque carte : un ressort entre le caisson (0) et la carte (i)
    for i in range(1, n):
        m_carte, amort_carte = cartes[i - 1][1], cartes[i - 1][5]
        k_carte = m_carte * (2 * np.pi * frequence_carte(cartes[i - 1])) ** 2
        c_carte = 2 * amort_carte * np.sqrt(k_carte * m_carte)
        K[0, 0] += k_carte
        K[i, i] += k_carte
        K[0, i] -= k_carte
        K[i, 0] -= k_carte
        C[0, 0] += c_carte
        C[i, i] += c_carte
        C[0, i] -= c_carte
        C[i, 0] -= c_carte

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
              f"Caisson : {masse_caisson} kg + {len(cartes)} cartes ; jeu disponible {jeu_mm} mm",
              "Critères : débattement en 4σ, flexion des cartes en 3σ au centre (Steinberg),",
              "           1re carte au moins 2 fois au-dessus des plots (règle de l'octave)", "",
              "Cartes (plaques FR4 en appui simple, fréquence calculée ; admissible de Steinberg) :"]
    for carte in cartes:
        lignes.append(f"    {carte[0]:20s} {carte[2]:.0f}x{carte[3]:.0f}x{carte[4]} mm, {carte[1]} kg"
                      f" -> {frequence_carte(carte):4.0f} Hz ; C = {carte[7]} -> admissible {admissible_steinberg_mm(carte):.3f} mm")
    lignes.append("")
    resultats = []  # pour chaque option : [T, grms, pire marge, critère limitant]

    for nom_plots, f_plots, amort_plots in options_plots:
        modes, z, T = reponse(f_plots, amort_plots, f)
        grms = [rms(abs(T[:, i]) ** 2 * dsp, f) for i in range(len(modes))]
        # z est en mètres PAR g d'entrée, la DSP en g²/Hz : |z|² x DSP x G² donne des m²/Hz
        debattement = 4 * rms(abs(z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # 4σ [mm]
        # règle de l'octave : la carte la plus basse doit être au moins 2 fois au-dessus des plots
        separation = min(frequence_carte(c) for c in cartes) / f_plots
        # MS = limite / calculé - 1 (pour la séparation, qui est un minimum : calculé / minimum - 1)
        marges_option = {"Séparation des modes": separation / 2.0 - 1,
                         "Débattement caisson": jeu_mm / debattement - 1}
        ms_sep, ms_deb = marges_option["Séparation des modes"], marges_option["Débattement caisson"]
        lignes += [f"--- {nom_plots} ({f_plots:.0f} Hz, amortissement {amort_plots:.0%})",
                   f"    Modes : {', '.join(f'{m:.0f}' for m in modes)} Hz",
                   f"    Séparation 1re carte / plots : x{separation:.2f} (minimum x2)  -> MS = {ms_sep:+.2f}  {verdict(ms_sep)}",
                   f"    Caisson : {grms[0]:.2f} Grms, débattement 4σ {debattement:.1f} mm  -> MS = {ms_deb:+.2f}  {verdict(ms_deb)}"]
        for i, carte in enumerate(cartes, start=1):
            # flexion de la carte par rapport au caisson, au centre de la carte, 3σ [mm]
            flexion = FACTEUR_CENTRE * 3 * rms(abs(z[:, i] - z[:, 0]) ** 2 * dsp * G**2, f) * 1000
            marges_option[carte[0]] = admissible_steinberg_mm(carte) / flexion - 1
            lignes.append(f"    {carte[0]:28s} {grms[i]:5.2f} Grms  flexion 3σ {flexion:.3f} mm  "
                          f"MS = {marges_option[carte[0]]:+6.2f}  {verdict(marges_option[carte[0]])}")
        limitant = min(marges_option, key=marges_option.get)  # le critère qui a la plus petite marge
        lignes += [f"    PIRE MARGE : {marges_option[limitant]:+.2f} ({limitant})", ""]
        resultats.append([T, grms, marges_option[limitant], limitant])

    lignes.append("CONCLUSION")
    for (nom_plots, _, _), (T, grms, pire, limitant) in zip(options_plots, resultats):
        lignes.append(f"  {nom_plots:28s} {verdict(pire):12s} limité par : {limitant:28s}"
                      f" pire carte = {max(grms[1:]) / grms_chassis:.1f} x châssis")
    conformes = [i for i in range(len(options_plots)) if resultats[i][2] >= 0]
    if conformes:
        meilleure_marge = max(conformes, key=lambda i: resultats[i][2])
        mieux_protegees = min(conformes, key=lambda i: max(resultats[i][1][1:]))
        lignes += [f"  -> Plus grande marge          : {options_plots[meilleure_marge][0]}",
                   f"  -> Cartes les mieux protégées : {options_plots[mieux_protegees][0]}"]
    else:
        lignes.append("  -> AUCUN choix de plots n'est conforme : il faut revoir le caisson ou les cartes")
    lignes += ["  Le choix final est un compromis : protéger les cartes (plots souples et amortis)",
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

    # figure de synthèse : 3 images = ce que font les plots, ce que reçoit chaque étage, le verdict
    noms = [o[0] for o in options_plots]
    figure, axes = plt.subplots(1, 3, figsize=(19, 5.8), gridspec_kw={"width_ratios": [1.4, 1.2, 1]})

    # 1. transmissibilité châssis -> caisson : les plots filtrent (< 1) ou amplifient (> 1) ?
    axe = axes[0]
    for nom, (T, grms, pire, limitant) in zip(noms, resultats):
        axe.loglog(f, abs(T[:, 0]), lw=2, label=nom)
    axe.axhline(1, color="grey", ls="--")
    axe.axvspan(profil_f[1], profil_f[2], color="tab:red", alpha=0.08)
    axe.text(np.sqrt(profil_f[1] * profil_f[2]), 12, "zone où la piste\nsecoue le plus", ha="center", fontsize=9, color="tab:red")
    axe.text(230, 8, "au-dessus de 1 : amplifie\nen dessous de 1 : isole", fontsize=9, color="grey")
    axe.set(title="Les plots amplifient-ils ou isolent-ils ?", xlabel="Fréquence [Hz]",
            ylabel="Mouvement du caisson / mouvement du châssis", xlim=(2, 1000), ylim=(5e-3, 30))
    axe.legend(fontsize=8, loc="lower left")
    axe.grid(True, which="both", alpha=0.3)

    # 2. niveau (Grms) reçu par chaque étage de la chaîne, pour chaque choix de plots
    etages = ["Châssis", "Caisson"] + [c[0].replace("Carte ", "Carte\n") for c in cartes]
    largeur = 0.8 / len(options_plots)
    for numero, (nom, (T, grms, pire, limitant)) in enumerate(zip(noms, resultats)):
        decalage = (numero - (len(options_plots) - 1) / 2) * largeur
        axes[1].bar(np.arange(len(etages)) + decalage, [grms_chassis] + grms, width=largeur, label=nom)
    axes[1].set_xticks(np.arange(len(etages)), etages, fontsize=9)
    axes[1].set(title="Niveau de vibration reçu à chaque étage", ylabel="Niveau [Grms]")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", alpha=0.3)

    # 3. verdict : la pire marge de chaque choix, et le critère qui limite
    pires = [r[2] for r in resultats]
    barres = axes[2].barh(noms, pires, color=["tab:green" if p >= 0 else "tab:red" for p in pires])
    for ligne, (pire, r) in enumerate(zip(pires, resultats)):
        axes[2].text(max(pire, 0) + 0.03, ligne, f"{pire:+.2f}\n{r[3]}", va="center", fontsize=8)
    axes[2].axvline(0, color="k")
    axes[2].set_xlim(min(min(pires), 0) - 0.2, max(pires) + 0.9)
    axes[2].set(title="Verdict : conforme si la marge > 0", xlabel="Pire marge de sécurité")
    axes[2].invert_yaxis()

    figure.suptitle("Radar sur camion : quels plots sous le caisson électronique ?  (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(os.path.join(dossier, "radar_camion.png"), dpi=150)
    plt.show()
