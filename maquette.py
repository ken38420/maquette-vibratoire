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
carte = [0.30, 200.0, 0.03]  # [masse kg, fréquence Hz, amortissement] de la carte dans le caisson

profil_f = [20, 80, 350, 2000]  # profil de vibration d'entrée : fréquences [Hz]
profil_dsp = [0.01, 0.04, 0.04, 0.007]  # et DSP [g²/Hz] (droites sur un graphe log-log)

jeu_mm = 5.0  # place disponible autour du caisson [mm]
bord, epaisseur, composant = 160 / 25.4, 1.6 / 25.4, 25 / 25.4  # carte, en pouces (pour Steinberg)

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
    debattement = 3 * rms(abs(z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # 3σ [mm]
    flexion = 3 * rms(abs(z[:, 1] - z[:, 0]) ** 2 * dsp * G**2, f) * 1000  # carte / caisson, 3σ [mm]
    admissible = 0.00022 * bord / (epaisseur * np.sqrt(composant)) * 25.4  # Steinberg [mm]
    separation = modes[1] / modes[0]  # règle de l'octave : au moins 2 (c'est un MINIMUM)
    return [["Séparation des modes f2/f1", separation, 2.0, separation / 2.0 - 1],
            ["Débattement silentblocs 3σ [mm]", debattement, jeu_mm, jeu_mm / debattement - 1],
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
    f_test = np.geomspace(0.1, 1e5, 100000)
    T_test = reponse([1.0, 100.0, 0.05], [1e-6, 1000.0, 0.05], f_test)[2]
    exact = np.sqrt(0.01 * np.pi * 100 * (1 + 4 * 0.05**2) / (4 * 0.05))
    ecart = abs(rms(abs(T_test[:, 0]) ** 2 * 0.01, f_test) / exact - 1)
    assert ecart < 1e-3, "le code ne retrouve pas la solution exacte"

    # note de calcul
    lignes = ["NOTE DE CALCUL - MAQUETTE VIBRATOIRE (valeurs illustratives)", "=" * 62, "",
              "1. HYPOTHÈSES",
              f"   Caisson : {caisson[0]} kg sur silentblocs à {caisson[1]} Hz, amortissement {caisson[2]:.0%}",
              f"   Carte   : {carte[0]} kg, fréquence propre {carte[1]} Hz, amortissement {carte[2]:.0%}",
              f"   Entrée  : {profil_f[0]}-{profil_f[-1]} Hz, {grms_entree:.2f} Grms, 1 axe vertical",
              "   Modèle  : 2 masses ponctuelles, ressorts et amortisseurs linéaires", "",
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

    # figure de synthèse
    figure, axes = plt.subplots(2, 2, figsize=(14, 9))
    axes[0, 0].loglog(f, abs(T[:, 0]), label="Caisson")
    axes[0, 0].loglog(f, abs(T[:, 1]), label="Carte")
    axes[0, 0].axhline(1, color="grey", ls="--")
    axes[0, 0].set(title="Transmissibilité (sortie / entrée)", xlabel="Fréquence [Hz]", xlim=(5, 2000))

    dans_profil = dsp > 0
    axes[0, 1].loglog(f[dans_profil], dsp[dans_profil], "k", lw=2, label=f"Entrée {grms_entree:.2f} Grms")
    axes[0, 1].loglog(f[dans_profil], (abs(T[:, 0]) ** 2 * dsp)[dans_profil], label=f"Caisson {grms_caisson:.2f} Grms")
    axes[0, 1].loglog(f[dans_profil], (abs(T[:, 1]) ** 2 * dsp)[dans_profil], label=f"Carte {grms_carte:.2f} Grms")
    axes[0, 1].set(title="DSP d'accélération [g²/Hz]", xlabel="Fréquence [Hz]")
    for axe in axes[0]:
        axe.legend()
        axe.grid(True, which="both", alpha=0.3)

    for axe, noms, valeurs, titre in [
        (axes[1, 0], [c[0] for c in criteres], [c[3] for c in criteres], "Marges, cas nominal (MS >= 0 : conforme)"),
        (axes[1, 1], [c[0] for c in cas], pires_ms, "Sensibilité : pire marge de chaque cas"),
    ]:
        axe.barh(noms, valeurs, color=["tab:green" if v >= 0 else "tab:red" for v in valeurs])
        axe.axvline(0, color="k")
        axe.set(title=titre, xlabel="MS")
    figure.suptitle("Maquette vibratoire : caisson sur silentblocs + carte (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(os.path.join(dossier, "synthese.png"), dpi=150)
    plt.show()
