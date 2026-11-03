"""Maquette vibratoire : caisson sur silentblocs + carte électronique.

Le caisson est posé sur des silentblocs, la carte est fixée dans le caisson.
Le porteur (camion, avion...) secoue la base avec une vibration aléatoire
définie par une DSP (densité spectrale de puissance, en g²/Hz).

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
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

type Vecteur = npt.NDArray[np.float64]

GRAVITE = 9.81  # [m/s²]
MM_PAR_POUCE = 25.4


class DonneeInvalide(Exception): ...


@dataclass
class Etage:
    """Un étage du système : une masse posée sur un ressort et un amortisseur."""

    masse_kg: float
    frequence_Hz: float  # fréquence propre de l'étage seul
    amortissement: float  # taux d'amortissement (0.10 = 10 %)


@dataclass
class Resultat:
    """Tout ce que le calcul produit pour un jeu de données."""

    frequences_propres_Hz: Vecteur
    deformees: npt.NDArray  # une colonne par mode, max = 1, caisson positif
    masses_effectives_kg: list[float]
    transmissibilite: npt.NDArray  # colonnes : caisson, carte
    dsp_sortie: npt.NDArray  # [g²/Hz], colonnes : caisson, carte
    grms_caisson: float
    grms_carte: float
    debattement_3sigma_mm: float  # caisson par rapport à la base
    deplacement_carte_3sigma_mm: float  # carte par rapport au caisson


# ============================================================
# 1. DONNÉES D'ENTRÉE (à modifier ici)
# ============================================================
CAISSON = Etage(masse_kg=25.0, frequence_Hz=25.0, amortissement=0.10)  # sur ses silentblocs
CARTE = Etage(masse_kg=0.30, frequence_Hz=200.0, amortissement=0.03)  # dans le caisson

# Profil d'entrée : points (fréquence [Hz], DSP [g²/Hz]). Entre deux points, la DSP suit une droite
# sur un graphe log-log, comme dans une spécification d'essai (+3 dB/octave, palier, -3 dB/octave).
PROFIL_FREQUENCES_HZ = [20, 80, 350, 2000]
PROFIL_DSP_G2_HZ = [0.01, 0.04, 0.04, 0.007]

JEU_MM = 5.0  # débattement disponible autour du caisson [mm]
CARTE_BORD_MM, CARTE_EPAISSEUR_MM, COMPOSANT_LONGUEUR_MM = 160.0, 1.6, 25.0  # pour Steinberg
RAPPORT_MODES_MIN = 2.0  # règle de l'octave : la carte vibre au moins 2 fois plus haut que la suspension

FREQUENCES_HZ = np.geomspace(1, 3000, 20000)  # fréquences de calcul


# ============================================================
# 2. CALCUL
# ============================================================
def valider(caisson: Etage, carte: Etage) -> None:
    """arrête le programme si une donnée n'a pas de sens physique"""
    for etage in (caisson, carte):
        if etage.masse_kg <= 0 or etage.frequence_Hz <= 0:
            raise DonneeInvalide(f"masse et fréquence doivent être > 0 : {etage}")
        if not 0 < etage.amortissement < 1:
            raise DonneeInvalide(f"l'amortissement doit être entre 0 et 1 : {etage}")


def calculer_dsp_entree(frequences: Vecteur) -> Vecteur:
    """DSP d'entrée [g²/Hz] : interpolation log-log du profil, nulle en dehors"""
    log_dsp = np.interp(np.log10(frequences), np.log10(PROFIL_FREQUENCES_HZ), np.log10(PROFIL_DSP_G2_HZ))
    dsp = 10**log_dsp
    hors_profil = (frequences < PROFIL_FREQUENCES_HZ[0]) | (frequences > PROFIL_FREQUENCES_HZ[-1])
    dsp[hors_profil] = 0.0
    return dsp


def calculer_rms(dsp: Vecteur, frequences: Vecteur) -> float:
    """valeur efficace (RMS) = racine de l'aire sous la DSP"""
    return float(np.sqrt(np.trapezoid(dsp, frequences)))


def construire_matrices(caisson: Etage, carte: Etage) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """matrices de masse, raideur et amortissement du système à 2 masses"""
    masse_suspendue = caisson.masse_kg + carte.masse_kg
    raideur_silentblocs = masse_suspendue * (2 * np.pi * caisson.frequence_Hz) ** 2  # k = m ω²
    raideur_carte = carte.masse_kg * (2 * np.pi * carte.frequence_Hz) ** 2
    amort_silentblocs = 2 * caisson.amortissement * np.sqrt(raideur_silentblocs * masse_suspendue)  # c = 2 ζ √(k m)
    amort_carte = 2 * carte.amortissement * np.sqrt(raideur_carte * carte.masse_kg)

    matrice_masse = np.array([[caisson.masse_kg, 0], [0, carte.masse_kg]])
    matrice_raideur = np.array(
        [[raideur_silentblocs + raideur_carte, -raideur_carte], [-raideur_carte, raideur_carte]]
    )
    matrice_amort = np.array([[amort_silentblocs + amort_carte, -amort_carte], [-amort_carte, amort_carte]])
    return matrice_masse, matrice_raideur, matrice_amort


def calculer_reponse(caisson: Etage, carte: Etage, frequences: Vecteur, dsp_entree: Vecteur) -> Resultat:
    """réponse du système quand la base est secouée par la DSP d'entrée"""
    valider(caisson, carte)
    matrice_masse, matrice_raideur, matrice_amort = construire_matrices(caisson, carte)
    un = np.ones(2)

    # modes propres : K φ = ω² M φ (les modes d'un système masse-ressort sont réels)
    omega_carre, deformees = np.linalg.eig(np.linalg.inv(matrice_masse) @ matrice_raideur)
    ordre = np.argsort(omega_carre.real)
    frequences_propres = np.sqrt(omega_carre.real[ordre]) / (2 * np.pi)
    deformees = deformees.real[:, ordre]
    # masse effective : part de la masse mise en mouvement par chaque mode
    masses_effectives = [float((mode @ matrice_masse @ un) ** 2 / (mode @ matrice_masse @ mode)) for mode in deformees.T]
    deformees = deformees / np.abs(deformees).max(axis=0) * np.sign(deformees[0])

    # réponse en fréquence : (K - ω²M + iωC) z = -M·1·a_base   (z = déplacement par rapport à la base)
    deplacements = np.array(
        [
            np.linalg.solve(matrice_raideur - omega**2 * matrice_masse + 1j * omega * matrice_amort, -matrice_masse @ un)
            for omega in 2 * np.pi * frequences
        ]
    )
    transmissibilite = 1 - (2 * np.pi * frequences[:, None]) ** 2 * deplacements  # accélération sortie / entrée

    dsp_sortie = np.abs(transmissibilite) ** 2 * dsp_entree[:, None]
    dsp_debattement = np.abs(deplacements[:, 0]) ** 2 * dsp_entree * GRAVITE**2  # [m²/Hz]
    dsp_carte = np.abs(deplacements[:, 1] - deplacements[:, 0]) ** 2 * dsp_entree * GRAVITE**2
    return Resultat(
        frequences_propres_Hz=frequences_propres,
        deformees=deformees,
        masses_effectives_kg=masses_effectives,
        transmissibilite=transmissibilite,
        dsp_sortie=dsp_sortie,
        grms_caisson=calculer_rms(dsp_sortie[:, 0], frequences),
        grms_carte=calculer_rms(dsp_sortie[:, 1], frequences),
        # 3σ : valeur dépassée seulement 0,3 % du temps
        debattement_3sigma_mm=3 * calculer_rms(dsp_debattement, frequences) * 1000,
        deplacement_carte_3sigma_mm=3 * calculer_rms(dsp_carte, frequences) * 1000,
    )


def calculer_deplacement_admissible_steinberg_mm() -> float:
    """critère de Steinberg (formule en pouces) : au-dessous, ~10 millions de cycles avant rupture des soudures"""
    bord, epaisseur, composant = (x / MM_PAR_POUCE for x in (CARTE_BORD_MM, CARTE_EPAISSEUR_MM, COMPOSANT_LONGUEUR_MM))
    return 0.00022 * bord / (epaisseur * np.sqrt(composant)) * MM_PAR_POUCE


def calculer_marges(resultat: Resultat) -> list[tuple[str, float, float, str, float]]:
    """critères : (nom, valeur calculée, limite, unité, marge de sécurité MS = admissible / calculé - 1)"""
    rapport_modes = resultat.frequences_propres_Hz[1] / resultat.frequences_propres_Hz[0]
    admissible_carte = calculer_deplacement_admissible_steinberg_mm()
    return [
        ("Séparation des modes f2/f1", rapport_modes, RAPPORT_MODES_MIN, "-", rapport_modes / RAPPORT_MODES_MIN - 1),
        ("Débattement silentblocs 3σ", resultat.debattement_3sigma_mm, JEU_MM, "mm",
         JEU_MM / resultat.debattement_3sigma_mm - 1),
        ("Déplacement carte 3σ (Steinberg)", resultat.deplacement_carte_3sigma_mm, admissible_carte, "mm",
         admissible_carte / resultat.deplacement_carte_3sigma_mm - 1),
    ]


def calculer_sensibilite(dsp_entree: Vecteur) -> list[tuple[str, Resultat, float]]:
    """refait le calcul quand les silentblocs changent ; renvoie (cas, résultat, pire marge)"""
    # la raideur change avec la température : fréquence x √(facteur de raideur)
    cas = [
        ("Nominal", 1.0, CAISSON.amortissement),
        ("Chaud : silentblocs -30 % raideur", 0.7, CAISSON.amortissement),
        ("Froid : silentblocs +50 % raideur", 1.5, CAISSON.amortissement),
        ("Amortissement faible (5 %)", 1.0, 0.05),
        ("Amortissement fort (20 %)", 1.0, 0.20),
    ]
    sensibilite = []
    for nom, facteur_raideur, amortissement in cas:
        caisson_modifie = Etage(CAISSON.masse_kg, CAISSON.frequence_Hz * np.sqrt(facteur_raideur), amortissement)
        resultat = calculer_reponse(caisson_modifie, CARTE, FREQUENCES_HZ, dsp_entree)
        pire_marge = min(critere[4] for critere in calculer_marges(resultat))
        sensibilite.append((nom, resultat, pire_marge))
    return sensibilite


def verifier_code() -> float:
    """contrôle sur un cas dont on connaît la réponse exacte :
    carte quasi sans masse -> le caisson est un oscillateur simple. Sous bruit blanc S0,
    Grms² = S0·π·fn·(1+4ζ²)/(4ζ)  (équation de Miles exacte)"""
    frequences = np.geomspace(0.1, 1e5, 100000)
    bruit_blanc = np.full_like(frequences, 0.01)
    resultat = calculer_reponse(Etage(1.0, 100.0, 0.05), Etage(1e-6, 1000.0, 0.05), frequences, bruit_blanc)
    grms_exact = np.sqrt(0.01 * np.pi * 100 * (1 + 4 * 0.05**2) / (4 * 0.05))
    ecart = abs(resultat.grms_caisson / grms_exact - 1)
    assert ecart < 1e-3, f"vérification échouée : écart {ecart:.1e}"
    return ecart


# ============================================================
# 3. SORTIES
# ============================================================
def conforme(marge: float) -> str:
    return "OK" if marge >= 0 else "NON CONFORME"


def ecrire_rapport(resultat, grms_entree, marges, sensibilite, ecart) -> str:
    """note de calcul au format texte"""
    masse_totale = CAISSON.masse_kg + CARTE.masse_kg
    amplification = np.abs(resultat.transmissibilite).max(axis=0)
    lignes = [
        "NOTE DE CALCUL - MAQUETTE VIBRATOIRE (valeurs illustratives)", "=" * 62, "",
        "1. HYPOTHÈSES",
        f"   Caisson : {CAISSON.masse_kg} kg sur silentblocs à {CAISSON.frequence_Hz} Hz, "
        f"amortissement {CAISSON.amortissement:.0%}",
        f"   Carte   : {CARTE.masse_kg} kg, fréquence propre {CARTE.frequence_Hz} Hz, "
        f"amortissement {CARTE.amortissement:.0%}",
        f"   Profil d'entrée : {PROFIL_FREQUENCES_HZ[0]}-{PROFIL_FREQUENCES_HZ[-1]} Hz, {grms_entree:.2f} Grms, 1 axe",
        "   Modèle : 2 masses ponctuelles, amortissement visqueux, comportement linéaire", "",
        "2. MODES PROPRES",
        "   Mode   Fréquence   Masse effective   Déformée caisson / carte",
    ]
    for numero in range(2):
        lignes.append(
            f"   {numero + 1}      {resultat.frequences_propres_Hz[numero]:7.1f} Hz   "
            f"{resultat.masses_effectives_kg[numero] / masse_totale:7.2%}            "
            f"{resultat.deformees[0, numero]:+.2f} / {resultat.deformees[1, numero]:+.2f}"
        )
    lignes += [
        "", "3. RÉPONSE ALÉATOIRE",
        f"   Entrée  : {grms_entree:.2f} Grms",
        f"   Caisson : {resultat.grms_caisson:.2f} Grms   (amplification max Q = {amplification[0]:.1f})",
        f"   Carte   : {resultat.grms_carte:.2f} Grms   (amplification max Q = {amplification[1]:.1f}, "
        f"pic 3σ = {3 * resultat.grms_carte:.1f} g)",
        "", "4. CRITÈRES  (marge de sécurité MS = admissible / calculé - 1, conforme si MS >= 0)",
    ]
    for nom, valeur, limite, unite, marge in marges:
        lignes.append(f"   {nom:34s} calculé {valeur:7.3f} {unite:2s}  limite {limite:7.3f} {unite:2s}  "
                      f"MS = {marge:+6.2f}  {conforme(marge)}")
    lignes += ["", "5. SENSIBILITÉ (pire marge parmi les 3 critères)"]
    for nom, resultat_cas, marge in sensibilite:
        lignes.append(f"   {nom:34s} f1 = {resultat_cas.frequences_propres_Hz[0]:5.1f} Hz  "
                      f"carte {resultat_cas.grms_carte:4.2f} Grms  débattement {resultat_cas.debattement_3sigma_mm:5.2f} mm  "
                      f"MS min = {marge:+6.2f}  {conforme(marge)}")
    pire_cas = min(sensibilite, key=lambda cas: cas[2])
    lignes += [
        "", "6. VÉRIFICATION DU CODE",
        f"   Oscillateur simple sous bruit blanc comparé à la solution exacte : écart {ecart:.1e}  OK", "",
        "CONCLUSION",
        f"   Pire cas : {pire_cas[0]} (MS = {pire_cas[2]:+.2f}) -> {conforme(pire_cas[2])}",
        "   Limites : 1 axe, masses ponctuelles, élastomère linéaire. À confirmer par un modèle",
        "   éléments finis puis par un essai au pot vibrant (voir recalage.py).",
    ]
    return "\n".join(lignes)


def tracer_synthese(resultat, dsp_entree, grms_entree, marges, sensibilite, chemin: str) -> None:
    """figure de synthèse en 4 graphiques"""
    frequences = FREQUENCES_HZ
    figure, axes = plt.subplots(2, 2, figsize=(14, 9))

    axe = axes[0, 0]
    axe.loglog(frequences, np.abs(resultat.transmissibilite[:, 0]), label="Caisson")
    axe.loglog(frequences, np.abs(resultat.transmissibilite[:, 1]), label="Carte")
    axe.axhline(1, color="grey", ls="--")
    axe.set(title="Transmissibilité (sortie / entrée)", xlabel="Fréquence [Hz]", xlim=(5, 2000))
    axe.legend()
    axe.grid(True, which="both", alpha=0.3)

    axe = axes[0, 1]
    dans_profil = dsp_entree > 0
    axe.loglog(frequences[dans_profil], dsp_entree[dans_profil], "k", lw=2, label=f"Entrée {grms_entree:.2f} Grms")
    axe.loglog(frequences[dans_profil], resultat.dsp_sortie[dans_profil, 0], label=f"Caisson {resultat.grms_caisson:.2f} Grms")
    axe.loglog(frequences[dans_profil], resultat.dsp_sortie[dans_profil, 1], label=f"Carte {resultat.grms_carte:.2f} Grms")
    axe.set(title="DSP d'accélération [g²/Hz]", xlabel="Fréquence [Hz]")
    axe.legend()
    axe.grid(True, which="both", alpha=0.3)

    for axe, noms, valeurs, titre in [
        (axes[1, 0], [m[0] for m in marges], [m[4] for m in marges], "Marges de sécurité, cas nominal (MS >= 0 : conforme)"),
        (axes[1, 1], [s[0] for s in sensibilite], [s[2] for s in sensibilite], "Sensibilité : pire marge de chaque cas"),
    ]:
        axe.barh(noms, valeurs, color=["tab:green" if valeur >= 0 else "tab:red" for valeur in valeurs])
        axe.axvline(0, color="k")
        axe.set(title=titre, xlabel="MS")

    figure.suptitle("Maquette vibratoire : caisson sur silentblocs + carte (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(chemin, dpi=150)


def main() -> None:
    dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
    os.makedirs(dossier, exist_ok=True)

    dsp_entree = calculer_dsp_entree(FREQUENCES_HZ)
    grms_entree = calculer_rms(dsp_entree, FREQUENCES_HZ)
    resultat = calculer_reponse(CAISSON, CARTE, FREQUENCES_HZ, dsp_entree)
    marges = calculer_marges(resultat)
    sensibilite = calculer_sensibilite(dsp_entree)
    ecart = verifier_code()

    rapport = ecrire_rapport(resultat, grms_entree, marges, sensibilite, ecart)
    print(rapport)
    with open(os.path.join(dossier, "rapport.txt"), "w", encoding="utf-8") as fichier:
        fichier.write(rapport)

    colonnes = [FREQUENCES_HZ, dsp_entree, resultat.dsp_sortie[:, 0], resultat.dsp_sortie[:, 1],
                np.abs(resultat.transmissibilite[:, 0]), np.abs(resultat.transmissibilite[:, 1])]
    np.savetxt(os.path.join(dossier, "courbes.csv"), np.column_stack(colonnes), delimiter=";", fmt="%.5g",
               comments="", header="f_Hz;DSP_entree_g2Hz;DSP_caisson_g2Hz;DSP_carte_g2Hz;T_caisson;T_carte")

    tracer_synthese(resultat, dsp_entree, grms_entree, marges, sensibilite, os.path.join(dossier, "synthese.png"))
    plt.show()


if __name__ == "__main__":
    main()
