"""
Maquette vibratoire : caisson sur silentblocs + carte électronique.

Le caisson (masse m1) est posé sur des silentblocs. La carte (masse m2) est fixée
dans le caisson. Le porteur (camion, avion...) secoue la base avec une vibration
aléatoire définie par une DSP (densité spectrale de puissance, en g²/Hz).

Le programme répond aux questions d'un bureau d'études :
  1. Quelles sont les fréquences propres ?                    -> modes
  2. La vibration est-elle amplifiée ou atténuée ?            -> transmissibilité
  3. Quel niveau reçoivent le caisson et la carte ?           -> Grms
  4. Les marges de sécurité sont-elles positives ?            -> critères
  5. Et si les silentblocs changent avec la température ?     -> sensibilité

Sorties dans le dossier "resultats" : rapport.txt, synthese.png, courbes.csv
Toutes les valeurs sont ILLUSTRATIVES (pas de données Thales).
"""
import os

import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# 1. DONNÉES D'ENTRÉE (à modifier ici)
# ============================================================
m1 = 25.0    # masse du caisson [kg]
f1 = 25.0    # fréquence de suspension sur silentblocs [Hz]
z1 = 0.10    # amortissement des silentblocs [-] (0.10 = 10 %)

m2 = 0.30    # masse de la carte électronique [kg]
f2 = 200.0   # fréquence propre de la carte [Hz]
z2 = 0.03    # amortissement de la carte [-]

# Profil de vibration d'entrée : points (fréquence [Hz], DSP [g²/Hz]).
# Entre deux points, la DSP suit une droite sur un graphe log-log,
# comme dans une spécification d'essai (ici +3 dB/octave, palier, -3 dB/octave).
profil_f = [20, 80, 350, 2000]
profil_S = [0.01, 0.04, 0.04, 0.007]

# Critères
jeu_mm = 5.0                # débattement disponible autour du caisson [mm]
B, h, L = 160.0, 1.6, 25.0  # carte : longueur du bord, épaisseur, longueur du plus gros composant [mm]
rapport_modes_min = 2.0     # règle de l'octave : la carte doit vibrer au moins 2 fois plus haut que la suspension

g = 9.81                            # [m/s²]
f = np.geomspace(1, 3000, 20000)    # fréquences de calcul [Hz]


# ============================================================
# 2. CALCUL
# ============================================================
def dsp_entree(f):
    """DSP d'entrée [g²/Hz] : interpolation log-log du profil, nulle en dehors."""
    S = 10 ** np.interp(np.log10(f), np.log10(profil_f), np.log10(profil_S))
    S[(f < profil_f[0]) | (f > profil_f[-1])] = 0
    return S


def rms(S, f):
    """Valeur efficace (RMS) = racine de l'aire sous la DSP."""
    return np.sqrt(np.trapezoid(S, f))


def calcul(f, S_in, m1, f1, z1, m2, f2, z2):
    """Réponse du système à 2 masses quand la base est secouée par la DSP S_in [g²/Hz]."""
    # Raideurs k = m (2πf)² et amortisseurs c = 2 ζ √(k m)
    k1 = (m1 + m2) * (2 * np.pi * f1) ** 2
    k2 = m2 * (2 * np.pi * f2) ** 2
    c1 = 2 * z1 * np.sqrt(k1 * (m1 + m2))
    c2 = 2 * z2 * np.sqrt(k2 * m2)
    M = np.array([[m1, 0], [0, m2]])
    K = np.array([[k1 + k2, -k2], [-k2, k2]])
    C = np.array([[c1 + c2, -c2], [-c2, c2]])
    un = np.ones(2)

    # Modes propres : K φ = ω² M φ
    w2, phi = np.linalg.eig(np.linalg.inv(M) @ K)
    w2, phi = w2.real, phi.real   # les modes d'un système masse-ressort sont réels
    ordre = np.argsort(w2)
    freqs = np.sqrt(w2[ordre]) / (2 * np.pi)
    phi = phi[:, ordre]
    # Masse effective : part de la masse totale mise en mouvement par chaque mode
    m_eff = [(p @ M @ un) ** 2 / (p @ M @ p) for p in phi.T]
    phi = phi / np.abs(phi).max(axis=0) * np.sign(phi[0])   # déformée : max = 1, caisson positif

    # Réponse en fréquence : (K - ω²M + iωC) z = -M·1·a_base   (z = déplacement par rapport à la base)
    z = np.array([np.linalg.solve(K - w**2 * M + 1j * w * C, -M @ un) for w in 2 * np.pi * f])
    T = 1 - (2 * np.pi * f[:, None]) ** 2 * z   # transmissibilité : accélération de sortie / entrée

    S_acc = np.abs(T) ** 2 * S_in[:, None]                    # DSP d'accélération [g²/Hz]
    S_debat = np.abs(z[:, 0]) ** 2 * S_in * g**2              # caisson par rapport à la base [m²/Hz]
    S_carte = np.abs(z[:, 1] - z[:, 0]) ** 2 * S_in * g**2    # carte par rapport au caisson [m²/Hz]
    return {
        "freqs": freqs, "phi": phi, "m_eff": m_eff, "T": T, "S_acc": S_acc,
        "grms_caisson": rms(S_acc[:, 0], f), "grms_carte": rms(S_acc[:, 1], f),
        "debat_3s_mm": 3 * rms(S_debat, f) * 1000,   # 3σ : valeur dépassée seulement 0,3 % du temps
        "carte_3s_mm": 3 * rms(S_carte, f) * 1000,
    }


def marges(r):
    """Tableau des critères : (nom, valeur calculée, limite, unité, marge de sécurité MS)."""
    po = 25.4  # mm par pouce : la formule de Steinberg est en pouces
    z_adm = 0.00022 * (B / po) / ((h / po) * np.sqrt(L / po)) * po
    rapport = r["freqs"][1] / r["freqs"][0]
    return [
        ("Séparation des modes f2/f1", rapport, rapport_modes_min, "-", rapport / rapport_modes_min - 1),
        ("Débattement silentblocs 3σ", r["debat_3s_mm"], jeu_mm, "mm", jeu_mm / r["debat_3s_mm"] - 1),
        ("Déplacement carte 3σ (Steinberg)", r["carte_3s_mm"], z_adm, "mm", z_adm / r["carte_3s_mm"] - 1),
    ]


def verifier():
    """Contrôle du code sur un cas dont on connaît la réponse exacte.
    Carte quasi sans masse -> le caisson est un oscillateur simple. Sous bruit blanc S0 :
    Grms² = S0·π·fn·(1+4ζ²)/(4ζ)  (équation de Miles exacte)."""
    fv = np.geomspace(0.1, 1e5, 100000)
    r = calcul(fv, np.full_like(fv, 0.01), 1.0, 100.0, 0.05, 1e-6, 1000.0, 0.05)
    exact = np.sqrt(0.01 * np.pi * 100 * (1 + 4 * 0.05**2) / (4 * 0.05))
    ecart = abs(r["grms_caisson"] / exact - 1)
    assert ecart < 1e-3, f"Vérification échouée : écart {ecart:.1e}"
    return ecart


S_in = dsp_entree(f)
grms_in = rms(S_in, f)
r = calcul(f, S_in, m1, f1, z1, m2, f2, z2)
crit = marges(r)
ecart = verifier()

# Sensibilité : la raideur des silentblocs change avec la température (fréquence x √facteur)
cas = [("Nominal", 1.0, z1),
       ("Chaud : silentblocs -30 % raideur", 0.7, z1),
       ("Froid : silentblocs +50 % raideur", 1.5, z1),
       ("Amortissement faible (5 %)", 1.0, 0.05),
       ("Amortissement fort (20 %)", 1.0, 0.20)]
sens = []
for nom, facteur, z in cas:
    rc = calcul(f, S_in, m1, f1 * np.sqrt(facteur), z, m2, f2, z2)
    sens.append((nom, rc, min(c[4] for c in marges(rc))))


# ============================================================
# 3. SORTIES
# ============================================================
dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
os.makedirs(dossier, exist_ok=True)

# --- Note de calcul (texte)
ok = lambda ms: "OK" if ms >= 0 else "NON CONFORME"
t = ["NOTE DE CALCUL - MAQUETTE VIBRATOIRE (valeurs illustratives)", "=" * 62, "",
     "1. HYPOTHÈSES",
     f"   Caisson : {m1} kg sur silentblocs à {f1} Hz, amortissement {z1:.0%}",
     f"   Carte   : {m2} kg, fréquence propre {f2} Hz, amortissement {z2:.0%}",
     f"   Profil d'entrée : {profil_f[0]}-{profil_f[-1]} Hz, {grms_in:.2f} Grms, 1 axe",
     "   Modèle : 2 masses ponctuelles, amortissement visqueux, comportement linéaire", "",
     "2. MODES PROPRES",
     "   Mode   Fréquence   Masse effective   Déformée caisson / carte"]
for j in range(2):
    t.append(f"   {j + 1}      {r['freqs'][j]:7.1f} Hz   {r['m_eff'][j] / (m1 + m2):7.2%}"
             f"            {r['phi'][0, j]:+.2f} / {r['phi'][1, j]:+.2f}")
T_max = np.abs(r["T"]).max(axis=0)
t += ["", "3. RÉPONSE ALÉATOIRE",
      f"   Entrée  : {grms_in:.2f} Grms",
      f"   Caisson : {r['grms_caisson']:.2f} Grms   (amplification max Q = {T_max[0]:.1f})",
      f"   Carte   : {r['grms_carte']:.2f} Grms   (amplification max Q = {T_max[1]:.1f}, pic 3σ = "
      f"{3 * r['grms_carte']:.1f} g)", "",
      "4. CRITÈRES  (marge de sécurité MS = admissible / calculé - 1, conforme si MS >= 0)"]
for nom, val, lim, u, ms in crit:
    t.append(f"   {nom:34s} calculé {val:7.3f} {u:2s}  limite {lim:7.3f} {u:2s}  MS = {ms:+6.2f}  {ok(ms)}")
t += ["", "5. SENSIBILITÉ (pire marge parmi les 3 critères)"]
for nom, rc, ms in sens:
    t.append(f"   {nom:34s} f1 = {rc['freqs'][0]:5.1f} Hz  carte {rc['grms_carte']:4.2f} Grms  "
             f"débattement {rc['debat_3s_mm']:5.2f} mm  MS min = {ms:+6.2f}  {ok(ms)}")
pire = min(sens, key=lambda s: s[2])
t += ["", "6. VÉRIFICATION DU CODE",
      f"   Oscillateur simple sous bruit blanc comparé à la solution exacte : écart {ecart:.1e}  OK", "",
      "CONCLUSION",
      f"   Pire cas : {pire[0]} (MS = {pire[2]:+.2f}) -> {ok(pire[2])}",
      "   Limites : 1 axe, masses ponctuelles, élastomère linéaire. À confirmer par un modèle",
      "   éléments finis puis par un essai au pot vibrant (recalage des fréquences et amortissements)."]
rapport = "\n".join(t)
print(rapport)
with open(os.path.join(dossier, "rapport.txt"), "w", encoding="utf-8") as fichier:
    fichier.write(rapport)

# --- Courbes pour Excel (séparateur ;)
np.savetxt(os.path.join(dossier, "courbes.csv"),
           np.column_stack([f, S_in, r["S_acc"][:, 0], r["S_acc"][:, 1], np.abs(r["T"][:, 0]), np.abs(r["T"][:, 1])]),
           delimiter=";", fmt="%.5g", comments="",
           header="f_Hz;DSP_entree_g2Hz;DSP_caisson_g2Hz;DSP_carte_g2Hz;T_caisson;T_carte")

# --- Figure de synthèse
fig, ax = plt.subplots(2, 2, figsize=(14, 9))
a = ax[0, 0]
a.loglog(f, np.abs(r["T"][:, 0]), label="Caisson")
a.loglog(f, np.abs(r["T"][:, 1]), label="Carte")
a.axhline(1, color="grey", ls="--")
a.set(title="Transmissibilité (sortie / entrée)", xlabel="Fréquence [Hz]", xlim=(5, 2000))
a.legend()
a.grid(True, which="both", alpha=0.3)

a = ax[0, 1]
bande = S_in > 0
a.loglog(f[bande], S_in[bande], "k", lw=2, label=f"Entrée {grms_in:.2f} Grms")
a.loglog(f[bande], r["S_acc"][bande, 0], label=f"Caisson {r['grms_caisson']:.2f} Grms")
a.loglog(f[bande], r["S_acc"][bande, 1], label=f"Carte {r['grms_carte']:.2f} Grms")
a.set(title="DSP d'accélération [g²/Hz]", xlabel="Fréquence [Hz]")
a.legend()
a.grid(True, which="both", alpha=0.3)

a = ax[1, 0]
ms = [c[4] for c in crit]
a.barh([c[0] for c in crit], ms, color=["tab:green" if m >= 0 else "tab:red" for m in ms])
a.axvline(0, color="k")
a.set(title="Marges de sécurité, cas nominal (MS >= 0 : conforme)", xlabel="MS")

a = ax[1, 1]
ms = [s[2] for s in sens]
a.barh([s[0] for s in sens], ms, color=["tab:green" if m >= 0 else "tab:red" for m in ms])
a.axvline(0, color="k")
a.set(title="Sensibilité : pire marge de chaque cas", xlabel="MS min")

fig.suptitle("Maquette vibratoire : caisson sur silentblocs + carte (valeurs illustratives)")
fig.tight_layout()
fig.savefig(os.path.join(dossier, "synthese.png"), dpi=150)
plt.show()
