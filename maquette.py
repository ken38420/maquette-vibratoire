"""Maquette vibratoire : un caisson sur silentblocs, avec une carte électronique dedans.

La base (le porteur) est secouée verticalement par une vibration aléatoire, donnée par
une DSP (densité spectrale de puissance, en g²/Hz).

Le programme calcule : les fréquences propres, la transmissibilité, les Grms,
les marges de sécurité, la sensibilité aux silentblocs, puis vérifie le code sur un cas exact.
Les résultats vont dans resultats/rapport.txt et resultats/synthese.png.
Valeurs ILLUSTRATIVES (ce ne sont pas des données Thales).
"""

import os

import matplotlib.pyplot as plt
import numpy as np

G = 9.81  # gravité [m/s²]

# ---------------- 1. DONNÉES (à modifier ici) ----------------
caisson = [25.0, 25.0, 0.10]  # [masse kg, fréquence Hz, amortissement] du caisson sur ses silentblocs
carte = [0.25, 200.0, 0.03]  # [masse kg, fréquence Hz, amortissement] de la carte dans le caisson
# (carte 160 x 100 x 2,0 mm de 0,25 kg : la formule de plaque de radar_camion.py donne 196 Hz, arrondi à 200)

profil_f = [20, 80, 350, 2000]  # profil de vibration d'entrée : fréquences [Hz]
profil_dsp = [0.01, 0.04, 0.04, 0.007]  # et DSP [g²/Hz] (droites sur un graphe log-log)

jeu_mm = 5.0  # place disponible autour du caisson [mm]
# carte pour Steinberg, en POUCES (la formule est en pouces) : longueur du bord, épaisseur, longueur du composant
bord, epaisseur, longueur_composant = 160 / 25.4, 2.0 / 25.4, 25 / 25.4

# La carte est réduite à UNE masse sur UN ressort : cela donne son déplacement "moyen" (modal).
# Au centre d'une plaque en appui sur ses 4 bords, la flexion est 16/π² ≈ 1,62 fois plus grande.
FACTEUR_CENTRE = 16 / np.pi**2

f = np.geomspace(1, 3000, 20000)  # fréquences de calcul [Hz]


# ---------------- 2. CALCUL ----------------
def dsp_profil(f):
    """DSP d'entrée [g²/Hz] : droites en log-log entre les points du profil, nulle en dehors"""
    dsp = 10 ** np.interp(np.log10(f), np.log10(profil_f), np.log10(profil_dsp))
    dsp[(f < profil_f[0]) | (f > profil_f[-1])] = 0
    return dsp


def rms(dsp, f):
    """valeur efficace (RMS) = racine de l'aire sous la DSP"""
    return np.sqrt(np.trapezoid(dsp, f))


def reponse(caisson, carte, f):
    """renvoie : fréquences propres [Hz], déplacements z et transmissibilité T (colonne 0 caisson, 1 carte)"""
    m1, f1, a1 = caisson
    m2, f2, a2 = carte
    k1 = (m1 + m2) * (2 * np.pi * f1) ** 2  # raideur des silentblocs : k = m ω²
    k2 = m2 * (2 * np.pi * f2) ** 2  # raideur de la carte
    c1 = 2 * a1 * np.sqrt(k1 * (m1 + m2))  # amortissement : c = 2 ζ √(k m)
    c2 = 2 * a2 * np.sqrt(k2 * m2)

    M = np.array([[m1, 0], [0, m2]])
    K = np.array([[k1 + k2, -k2], [-k2, k2]])
    C = np.array([[c1 + c2, -c2], [-c2, c2]])

    modes = np.sort(np.sqrt(np.linalg.eigvals(np.linalg.inv(M) @ K).real)) / (2 * np.pi)
    w = 2 * np.pi * f
    # à chaque fréquence : (K - ω²M + iωC) z = -M·1·a   (z = déplacement par rapport à la base)
    z = np.array([np.linalg.solve(K - wi**2 * M + 1j * wi * C, -M @ [1, 1]) for wi in w])
    T = 1 - w[:, None] ** 2 * z  # accélération de la pièce / accélération de la base
    return modes, z, T


def marges(caisson):
    """les 3 critères : [nom, valeur calculée, limite, marge de sécurité MS] (conforme si MS >= 0)"""
    modes, z, T = reponse(caisson, carte, f)
    dsp = dsp_profil(f)
    # z est en mètres PAR g d'entrée, et la DSP est en g²/Hz : |z|² x DSP x G² donne des m²/Hz
    # débattement : 4σ, car sur une heure d'essai le plus grand pic atteint environ 4σ (loi de Rayleigh)
    debattement = 4 * rms(abs(z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # [mm]
    # flexion de la carte par rapport au caisson, au centre de la carte, 3σ (convention de Steinberg)
    flexion = FACTEUR_CENTRE * 3 * rms(abs(z[:, 1] - z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # [mm]
    # Steinberg : Z = 0,00022 B / (C h r √L) en pouces, pour environ 20 millions d'alternances.
    # Ici C = 1 (composant à broches) et r = 1 (composant au centre, le cas le plus sévère).
    admissible = 0.00022 * bord / (epaisseur * np.sqrt(longueur_composant)) * 25.4  # [mm]
    separation = modes[1] / modes[0]  # règle de l'octave : au moins 2 (c'est un MINIMUM)
    return [["Séparation des modes f2/f1", separation, 2.0, separation / 2.0 - 1],
            ["Débattement silentblocs 4σ [mm]", debattement, jeu_mm, jeu_mm / debattement - 1],
            ["Flexion carte 3σ [mm] (Steinberg)", flexion, admissible, admissible / flexion - 1]]


def verdict(ms):
    return "OK" if ms >= 0 else "NON CONFORME"


# ---------------- 3. PROGRAMME PRINCIPAL ----------------
if __name__ == "__main__":
    dsp = dsp_profil(f)
    modes, z, T = reponse(caisson, carte, f)
    grms_entree = rms(dsp, f)
    grms_caisson = rms(abs(T[:, 0]) ** 2 * dsp, f)
    grms_carte = rms(abs(T[:, 1]) ** 2 * dsp, f)
    criteres = marges(caisson)

    # sensibilité : les silentblocs changent avec la température (fréquence x √facteur de raideur)
    m, fr, a = caisson
    cas = [["Nominal", caisson],
           ["Chaud : silentblocs -30 % raideur", [m, fr * np.sqrt(0.7), a]],
           ["Froid : silentblocs +50 % raideur", [m, fr * np.sqrt(1.5), a]],
           ["Amortissement faible (5 %)", [m, fr, 0.05]],
           ["Amortissement fort (20 %)", [m, fr, 0.20]]]
    pires_ms = [min(c[3] for c in marges(caisson_cas)) for nom, caisson_cas in cas]

    # vérification du code : carte presque sans masse -> oscillateur simple sous bruit blanc,
    # dont la réponse exacte est connue : Grms² = S0·π·fn·(1+4ζ²)/(4ζ)
    # (Miles est l'approximation de cette formule pour ζ petit : Grms² ≈ π/2 · fn · Q · S0, avec Q = 1/2ζ)
    S0, fn, zeta = 0.01, 100.0, 0.05  # bruit blanc [g²/Hz], fréquence [Hz], amortissement
    f_test = np.geomspace(0.1, 1e5, 100000)
    T_test = reponse([1.0, fn, zeta], [1e-6, 1000.0, zeta], f_test)[2]
    exact = np.sqrt(S0 * np.pi * fn * (1 + 4 * zeta**2) / (4 * zeta))
    ecart = abs(rms(abs(T_test[:, 0]) ** 2 * S0, f_test) / exact - 1)
    assert ecart < 1e-3, "le code ne retrouve pas la solution exacte"

    # note de calcul
    lignes = ["NOTE DE CALCUL - MAQUETTE VIBRATOIRE (valeurs illustratives)", "=" * 62, "",
              "1. HYPOTHÈSES",
              f"   Caisson : {caisson[0]} kg sur silentblocs à {caisson[1]} Hz, amortissement {caisson[2]:.0%}",
              f"   Carte   : {carte[0]} kg, fréquence propre {carte[1]} Hz, amortissement {carte[2]:.0%}",
              f"   Entrée  : {profil_f[0]}-{profil_f[-1]} Hz, {grms_entree:.2f} Grms, 1 axe vertical",
              "   Modèle  : 2 masses ponctuelles, ressorts et amortisseurs linéaires",
              "   Critères: débattement en 4σ (plus grand pic sur 1 h d'essai), flexion de la carte",
              "             en 3σ au centre de la carte (Steinberg, C = 1, r = 1)", "",
              f"2. FRÉQUENCES PROPRES : {modes[0]:.1f} Hz et {modes[1]:.1f} Hz", "",
              "3. NIVEAUX REÇUS",
              f"   Caisson : {grms_caisson:.2f} Grms (amplification max x{abs(T[:, 0]).max():.1f})",
              f"   Carte   : {grms_carte:.2f} Grms (amplification max x{abs(T[:, 1]).max():.1f})", "",
              "4. CRITÈRES (marge de sécurité MS, conforme si MS >= 0)"]
    for nom, valeur, limite, ms in criteres:
        lignes.append(f"   {nom:36s} calculé {valeur:7.3f}  limite {limite:7.3f}  MS = {ms:+6.2f}  {verdict(ms)}")
    lignes += ["", "5. SENSIBILITÉ (pire marge des 3 critères)"]
    for (nom, caisson_cas), ms in zip(cas, pires_ms):
        lignes.append(f"   {nom:36s} f1 = {caisson_cas[1]:5.1f} Hz  MS min = {ms:+6.2f}  {verdict(ms)}")
    lignes += ["", "6. VÉRIFICATION DU CODE",
               f"   Oscillateur simple sous bruit blanc / solution exacte : écart {ecart:.1e}  OK", "",
               "CONCLUSION",
               f"   Pire cas : {cas[np.argmin(pires_ms)][0]} (MS = {min(pires_ms):+.2f}) -> {verdict(min(pires_ms))}",
               "   Limites : 1 axe, masses ponctuelles, élastomère linéaire. À confirmer par un modèle",
               "   éléments finis puis par un essai au pot vibrant (voir recalage.py)."]
    rapport = "\n".join(lignes)
    print(rapport)
    dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
    os.makedirs(dossier, exist_ok=True)
    with open(os.path.join(dossier, "rapport.txt"), "w", encoding="utf-8") as fichier:
        fichier.write(rapport)

    # figure de synthèse : 3 images, de gauche à droite = ce qui se passe, est-ce conforme, est-ce robuste
    figure, axes = plt.subplots(1, 3, figsize=(19, 5.5), gridspec_kw={"width_ratios": [1.5, 1, 1]})

    # 1. transmissibilité : combien chaque pièce bouge par rapport à la base, selon la fréquence
    axe = axes[0]
    axe.loglog(f, abs(T[:, 0]), lw=2, label="Caisson")
    axe.loglog(f, abs(T[:, 1]), lw=2, label="Carte")
    axe.axhline(1, color="grey", ls="--")
    axe.annotate(f"résonance des silentblocs\n({modes[0]:.0f} Hz, x{abs(T[:, 0]).max():.0f})",
                 xy=(modes[0], abs(T[:, 0]).max()), xytext=(45, 9), fontsize=9,
                 arrowprops=dict(arrowstyle="->", color="k"))
    axe.annotate(f"résonance de la carte\n({modes[1]:.0f} Hz)", xy=(modes[1], abs(T[:, 1])[np.argmin(abs(f - modes[1]))]),
                 xytext=(320, 2), fontsize=9, arrowprops=dict(arrowstyle="->", color="k"))
    axe.text(6, 0.004, "au-dessus de 1 : amplifie\nen dessous de 1 : isole", fontsize=9, color="grey")
    axe.set(title="Ce qui bouge, et à quelle fréquence", xlabel="Fréquence [Hz]",
            ylabel="Mouvement de la pièce / mouvement de la base", xlim=(5, 2000), ylim=(1e-3, 20))
    axe.legend(loc="lower left")
    axe.grid(True, which="both", alpha=0.3)

    # 2. les 3 critères du cas nominal (vert = conforme, rouge = non conforme)
    ms_nominal = [c[3] for c in criteres]
    noms_courts = ["Modes assez séparés ?", "Le caisson ne tape pas ?", "La carte tient (Steinberg) ?"]
    barres = axes[1].barh(noms_courts, ms_nominal, color=["tab:green" if v >= 0 else "tab:red" for v in ms_nominal])
    axes[1].bar_label(barres, fmt="%+.2f", padding=3)
    axes[1].axvline(0, color="k")
    axes[1].set(title="Les 3 critères (marge > 0 : conforme)", xlabel="Marge de sécurité",
                xlim=(min(ms_nominal + [0]) - 0.2, max(ms_nominal) * 1.25))
    axes[1].invert_yaxis()

    # 3. robustesse : la pire marge reste-t-elle positive si les silentblocs changent ?
    barres = axes[2].barh([c[0] for c in cas], pires_ms,
                          color=["tab:green" if v >= 0 else "tab:red" for v in pires_ms])
    axes[2].bar_label(barres, fmt="%+.2f", padding=3)
    axes[2].axvline(0, color="k")
    axes[2].set(title="Et si les silentblocs changent ? (pire marge)", xlabel="Pire marge de sécurité",
                xlim=(min(pires_ms + [0]) - 0.2, max(pires_ms) * 1.25))
    axes[2].invert_yaxis()

    figure.suptitle(f"Caisson de {caisson[0]:.0f} kg sur silentblocs à {caisson[1]:.0f} Hz, carte à {carte[1]:.0f} Hz, "
                    f"entrée {grms_entree:.1f} Grms  (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(os.path.join(dossier, "synthese.png"), dpi=150)
    plt.show()
