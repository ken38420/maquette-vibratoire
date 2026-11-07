"""Recalage du modèle sur un essai au pot vibrant, par descente de gradient.

La fiche technique des silentblocs dit 25 Hz et 10 % d'amortissement, mais le vrai produit
ne se comporte pas exactement comme ça. On mesure donc la transmissibilité du caisson
(balayage sinus à bas niveau), puis on cherche les paramètres du modèle qui collent le mieux
à la mesure : c'est le RECALAGE.

Méthode du cours (réseaux de neurones > perte et entraînement) :
  - perte J = erreur quadratique moyenne entre modèle et mesure (avec le facteur 1/2)
  - descente de gradient : θ <- θ - α ∇J, arrêt quand le pas devient minuscule
  - 70 % des mesures servent à recaler, 30 % à vérifier (jeu de test)

L'essai est SIMULÉ : on part d'un "vrai" caisson caché et on ajoute 5 % de bruit.
Limite honnête : les points de test viennent du même balayage que ceux d'entraînement.
Une vraie validation se ferait sur un AUTRE essai (autre niveau, autre axe, autre exemplaire).
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from maquette import caisson, carte, marges, reponse

vrai_caisson = [25.0, 22.0, 0.14]  # inconnu en vrai : c'est ce qu'on cherche
f_essai = np.geomspace(10, 150, 60)  # points de mesure [Hz]
alpha = 0.5  # taille du pas
rng = np.random.default_rng(0)  # graine fixe : le même bruit à chaque lancement


def vers_caisson(p):
    """paramètres θ = [fréquence / 25 Hz, amortissement / 0.10] -> caisson"""
    return [caisson[0], p[0] * 25.0, p[1] * 0.10]


def predire(p, f):
    """transmissibilité du caisson prédite par le modèle"""
    return abs(reponse(vers_caisson(p), carte, f)[2][:, 0])


def perte(p, f, mesures):
    """J = moyenne de 1/2 (prédiction - mesure)², en log car T varie sur plusieurs décades"""
    return np.mean(0.5 * (np.log10(predire(p, f)) - np.log10(mesures)) ** 2)


def gradient(p, f, mesures):
    """∇J par différences finies : (J(θ+h) - J(θ-h)) / 2h pour chaque paramètre"""
    h = 1e-6
    return np.array([(perte(p + h * e, f, mesures) - perte(p - h * e, f, mesures)) / (2 * h) for e in np.eye(2)])


# mesure simulée, puis séparation entraînement / test
mesures = abs(reponse(vrai_caisson, carte, f_essai)[2][:, 0]) * (1 + 0.05 * rng.standard_normal(f_essai.size))
melange = rng.permutation(f_essai.size)
entr, test = np.sort(melange[:42]), np.sort(melange[42:])

# descente de gradient, en partant de la fiche technique
p_initial = np.array([caisson[1] / 25.0, caisson[2] / 0.10])
p = p_initial.copy()
historique = [perte(p, f_essai[entr], mesures[entr])]
for iteration in range(2000):
    pas = alpha * gradient(p, f_essai[entr], mesures[entr])
    p = p - pas
    historique.append(perte(p, f_essai[entr], mesures[entr]))
    if np.linalg.norm(pas) < 1e-7:  # critère d'arrêt : le pas est devenu minuscule
        break

initial, recale = vers_caisson(p_initial), vers_caisson(p)
print("RECALAGE DU MODÈLE SUR L'ESSAI (essai simulé)")
print(f"  Itérations : {len(historique) - 1}")
print("                  Fiche technique   Recalé   Vrai (caché)")
print(f"  Fréquence [Hz]  {initial[1]:10.2f}   {recale[1]:10.2f}   {vrai_caisson[1]:8.2f}")
print(f"  Amortissement   {initial[2]:10.1%}   {recale[2]:10.1%}   {vrai_caisson[2]:8.1%}")
print(f"  Perte J sur le jeu de TEST : {perte(p_initial, f_essai[test], mesures[test]):.2e}"
      f" -> {perte(p, f_essai[test], mesures[test]):.2e}")
print(f"  Pire marge de sécurité     : {min(c[3] for c in marges(initial)):+.2f}"
      f" -> {min(c[3] for c in marges(recale)):+.2f}")

figure, axes = plt.subplots(1, 2, figsize=(13, 5))
f_fin = np.geomspace(10, 150, 400)
axes[0].loglog(f_essai[entr], mesures[entr], "o", label="Mesures (entraînement)")
axes[0].loglog(f_essai[test], mesures[test], "s", mfc="none", label="Mesures (test)")
axes[0].loglog(f_fin, predire(p_initial, f_fin), "--", label="Modèle fiche technique")
axes[0].loglog(f_fin, predire(p, f_fin), "k", label="Modèle recalé")
axes[0].set(title="Transmissibilité du caisson : essai et modèle", xlabel="Fréquence [Hz]", ylabel="|T| [-]")
axes[0].legend()
axes[0].grid(True, which="both", alpha=0.3)
axes[1].semilogy(historique)
axes[1].set(title="Descente de gradient : perte J à chaque itération", xlabel="Itération", ylabel="J")
axes[1].grid(True, which="both", alpha=0.3)
figure.tight_layout()
dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
os.makedirs(dossier, exist_ok=True)  # crée le dossier s'il n'existe pas encore
figure.savefig(os.path.join(dossier, "recalage.png"), dpi=150)
plt.show()
