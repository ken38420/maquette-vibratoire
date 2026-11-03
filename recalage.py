"""Recalage du modèle sur un essai au pot vibrant, par descente de gradient.

Situation réelle : la fiche technique des silentblocs donne 25 Hz et 10 % d'amortissement,
mais le vrai produit ne se comporte pas exactement comme ça. On fait donc un essai
(balayage sinus à bas niveau) qui mesure la transmissibilité du caisson, puis on cherche
les paramètres du modèle qui reproduisent le mieux la mesure : c'est le RECALAGE.

Méthode, reprise du cours (Forward Neural Network > Loss / Training) :
  - perte J = erreur quadratique moyenne entre modèle et mesure (MSE avec le facteur 1/2)
  - descente de gradient : θ(k+1) = θ(k) - α ∇J(θ(k)), arrêt quand ||α ∇J|| < ε
  - les mesures sont séparées en jeu d'entraînement (pour recaler) et jeu de test (pour juger)
Ici θ ne contient que 2 paramètres (fréquence et amortissement des silentblocs) :
le gradient est calculé par différences finies, l'autodifférentiation du cours n'est pas utile.

L'essai est SIMULÉ (pas de vraie mesure) : on part d'un "vrai" caisson caché et on ajoute du bruit.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from maquette import CAISSON, CARTE, FREQUENCES_HZ, Etage, Vecteur, calculer_dsp_entree, calculer_marges, calculer_reponse

SEED = 0  # graine du hasard : le même bruit à chaque lancement (répétabilité)
VRAI_CAISSON = Etage(masse_kg=25.0, frequence_Hz=22.0, amortissement=0.14)  # inconnu en vrai : c'est ce qu'on cherche
BRUIT_MESURE = 0.05  # 5 % de bruit sur chaque point de mesure
FREQUENCES_ESSAI_HZ = np.geomspace(10, 150, 60)  # points du balayage sinus
PART_ENTRAINEMENT = 0.7  # 70 % des points pour recaler, 30 % pour tester

# Paramètres θ normalisés pour qu'ils soient du même ordre de grandeur (~1)
FREQUENCE_REFERENCE_HZ = 25.0
AMORTISSEMENT_REFERENCE = 0.10

LEARNING_RATE = 0.5  # α : taille du pas
NB_ITERATIONS_MAX = 2000
EPSILON = 1e-7  # critère d'arrêt sur ||α ∇J||


def construire_caisson(parametres: Vecteur) -> Etage:
    """θ = [fréquence / 25 Hz, amortissement / 0.10] -> caisson"""
    return Etage(CAISSON.masse_kg, parametres[0] * FREQUENCE_REFERENCE_HZ, parametres[1] * AMORTISSEMENT_REFERENCE)


def predire(parametres: Vecteur, frequences: Vecteur) -> Vecteur:
    """transmissibilité du caisson prédite par le modèle"""
    resultat = calculer_reponse(construire_caisson(parametres), CARTE, frequences, np.zeros_like(frequences))
    return np.abs(resultat.transmissibilite[:, 0])


def simuler_essai(rng: np.random.Generator) -> Vecteur:
    """mesure simulée : le vrai caisson + du bruit de mesure"""
    resultat = calculer_reponse(VRAI_CAISSON, CARTE, FREQUENCES_ESSAI_HZ, np.zeros_like(FREQUENCES_ESSAI_HZ))
    transmissibilite_vraie = np.abs(resultat.transmissibilite[:, 0])
    return transmissibilite_vraie * (1 + BRUIT_MESURE * rng.standard_normal(transmissibilite_vraie.size))


def calculer_perte(parametres: Vecteur, frequences: Vecteur, mesures: Vecteur) -> float:
    """J = 1/n Σ 1/2 (prédiction - mesure)², en log car la transmissibilité varie sur plusieurs décades"""
    ecarts = np.log10(predire(parametres, frequences)) - np.log10(mesures)
    return float(np.mean(0.5 * ecarts**2))


def calculer_gradient(parametres: Vecteur, frequences: Vecteur, mesures: Vecteur) -> Vecteur:
    """∇J par différences finies centrées : (J(θ+h) - J(θ-h)) / 2h pour chaque paramètre"""
    pas = 1e-6
    gradient = np.zeros_like(parametres)
    for indice in range(parametres.size):
        decalage = np.zeros_like(parametres)
        decalage[indice] = pas
        perte_plus = calculer_perte(parametres + decalage, frequences, mesures)
        perte_moins = calculer_perte(parametres - decalage, frequences, mesures)
        gradient[indice] = (perte_plus - perte_moins) / (2 * pas)
    return gradient


def recaler(parametres_initiaux: Vecteur, frequences: Vecteur, mesures: Vecteur) -> tuple[Vecteur, list[float]]:
    """descente de gradient ; renvoie les paramètres recalés et l'historique de la perte"""
    parametres = parametres_initiaux.copy()
    historique = [calculer_perte(parametres, frequences, mesures)]
    for _ in range(NB_ITERATIONS_MAX):
        avance = LEARNING_RATE * calculer_gradient(parametres, frequences, mesures)
        parametres = parametres - avance
        historique.append(calculer_perte(parametres, frequences, mesures))
        if np.linalg.norm(avance) < EPSILON:
            break
    return parametres, historique


def main() -> None:
    rng = np.random.default_rng(SEED)
    mesures = simuler_essai(rng)

    # séparation entraînement / test (tirage au hasard des points de mesure)
    melange = rng.permutation(FREQUENCES_ESSAI_HZ.size)
    nb_entrainement = int(PART_ENTRAINEMENT * FREQUENCES_ESSAI_HZ.size)
    entrainement, test = np.sort(melange[:nb_entrainement]), np.sort(melange[nb_entrainement:])

    parametres_initiaux = np.array([CAISSON.frequence_Hz / FREQUENCE_REFERENCE_HZ,
                                    CAISSON.amortissement / AMORTISSEMENT_REFERENCE])  # fiche technique
    parametres, historique = recaler(parametres_initiaux, FREQUENCES_ESSAI_HZ[entrainement], mesures[entrainement])
    caisson_initial, caisson_recale = construire_caisson(parametres_initiaux), construire_caisson(parametres)

    perte_test_initiale = calculer_perte(parametres_initiaux, FREQUENCES_ESSAI_HZ[test], mesures[test])
    perte_test_finale = calculer_perte(parametres, FREQUENCES_ESSAI_HZ[test], mesures[test])

    # les marges changent-elles avec le modèle recalé ?
    dsp_entree = calculer_dsp_entree(FREQUENCES_HZ)
    pire_marge = {}
    for nom, caisson in [("initial", caisson_initial), ("recalé", caisson_recale)]:
        resultat = calculer_reponse(caisson, CARTE, FREQUENCES_HZ, dsp_entree)
        pire_marge[nom] = min(critere[4] for critere in calculer_marges(resultat))

    print("RECALAGE DU MODÈLE SUR L'ESSAI (essai simulé)")
    print(f"  Itérations de descente de gradient : {len(historique) - 1}")
    print("                     Fiche technique   Recalé    Vrai (caché)")
    print(f"  Fréquence [Hz]     {caisson_initial.frequence_Hz:10.2f}   {caisson_recale.frequence_Hz:10.2f}"
          f"   {VRAI_CAISSON.frequence_Hz:10.2f}")
    print(f"  Amortissement      {caisson_initial.amortissement:10.1%}   {caisson_recale.amortissement:10.1%}"
          f"   {VRAI_CAISSON.amortissement:10.1%}")
    print(f"  Perte J sur le jeu de TEST : {perte_test_initiale:.2e} -> {perte_test_finale:.2e}")
    print(f"  Pire marge de sécurité     : {pire_marge['initial']:+.2f} -> {pire_marge['recalé']:+.2f}")

    figure, axes = plt.subplots(1, 2, figsize=(13, 5))
    frequences_fines = np.geomspace(10, 150, 400)
    axe = axes[0]
    axe.loglog(FREQUENCES_ESSAI_HZ[entrainement], mesures[entrainement], "o", label="Mesures (entraînement)")
    axe.loglog(FREQUENCES_ESSAI_HZ[test], mesures[test], "s", mfc="none", label="Mesures (test)")
    axe.loglog(frequences_fines, predire(parametres_initiaux, frequences_fines), "--", label="Modèle fiche technique")
    axe.loglog(frequences_fines, predire(parametres, frequences_fines), "k", label="Modèle recalé")
    axe.set(title="Transmissibilité du caisson : essai et modèle", xlabel="Fréquence [Hz]")
    axe.legend()
    axe.grid(True, which="both", alpha=0.3)

    axe = axes[1]
    axe.semilogy(historique)
    axe.set(title="Descente de gradient : perte J à chaque itération", xlabel="Itération", ylabel="J (entraînement)")
    axe.grid(True, which="both", alpha=0.3)

    figure.tight_layout()
    dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
    os.makedirs(dossier, exist_ok=True)
    figure.savefig(os.path.join(dossier, "recalage.png"), dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
