"""Matrice d'essais : mutualiser la qualification environnementale de plusieurs projets.

C'est la mission "standardiser le référentiel" de l'offre :
  1. on recense les exigences de chaque projet (une ligne par projet) ;
  2. on prend l'ENVELOPPE (le cas le plus sévère de chaque environnement) ;
  3. on regarde si un produit DÉJÀ QUALIFIÉ couvre déjà cette enveloppe : sa preuve est réutilisable
     (qualification par similitude, pas de nouvel essai) ;
  4. on compte les essais : un par projet (sans mutualisation) ou un seul essai enveloppe (avec).

Les résultats vont dans resultats/matrice_essais.txt et resultats/matrice_essais.png.
Valeurs ILLUSTRATIVES (ce ne sont pas des données Thales). Les numéros de méthode AECTP sont
ceux de l'édition que je connais : à vérifier dans l'édition en vigueur.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

# ---------------- 1. DONNÉES (à modifier ici) ----------------
# environnements : [nom, norme / méthode, unité, sens]
# sens "max" : la valeur la plus haute est la plus sévère ; "min" : la plus basse ; "oui" : exigé ou pas
environnements = [["Vibration aléatoire", "AECTP-400 m. 401", "Grms", "max"],
                  ["Choc", "AECTP-400 m. 403", "g", "max"],
                  ["Chaleur (en marche)", "AECTP-300 m. 302", "°C", "max"],
                  ["Froid (en marche)", "AECTP-300 m. 303", "°C", "min"],
                  ["Chaleur humide", "AECTP-300 m. 306", "jours", "max"],
                  ["Brouillard salin", "AECTP-300 m. 309", "h", "max"],
                  ["CEM rayonnée", "AECTP-500 / 461 RS103", "V/m", "max"],
                  ["Feu", "hors AECTP", "", "oui"]]

# exigences de chaque projet, dans l'ordre des environnements (None = pas d'exigence)
projets = [["Radar fixe en shelter", [1.0, 15, 55, -33, 10, 48, 50, True]],
           ["Radar sur camion", [1.6, 30, 55, -40, 10, 96, 50, None]],
           ["Version aérotransportable", [2.5, 20, 49, -46, None, None, 200, None]]]

# produit déjà qualifié : niveaux qu'il a déjà tenus en essai (sa preuve peut resservir).
# La similitude se justifie : même conception, mêmes matériaux, même fixation. Sinon, il faut un essai.
deja_qualifie = ["Caisson de la génération précédente", [2.0, 30, 55, -40, 10, 48, 100, None]]

# ordre conseillé des essais sur UNE seule éprouvette : mécanique d'abord (il révèle les
# fissures et desserrages), le climatique ensuite, le brouillard salin en dernier (il corrode).
# Le feu n'y est pas : c'est un essai destructif, fait sur une éprouvette (ou un matériau) à part.
ordre_essais = [0, 1, 6, 3, 2, 4, 5]


# ---------------- 2. CALCUL ----------------
def plus_severe(valeurs, sens):
    """enveloppe d'un environnement : la valeur la plus sévère parmi celles exigées"""
    exigees = [v for v in valeurs if v is not None]
    if not exigees:
        return None
    if sens == "max":
        return max(exigees)
    if sens == "min":
        return min(exigees)
    return True  # "oui" : au moins un projet l'exige


def couvert(preuve, besoin, sens):
    """le niveau déjà prouvé couvre-t-il le besoin ?"""
    if besoin is None:
        return True
    if preuve is None:
        return False
    if sens == "max":
        return preuve >= besoin
    if sens == "min":
        return preuve <= besoin
    return bool(preuve)


def texte(valeur):
    """valeur affichée dans la matrice"""
    if valeur is None:
        return "-"
    if valeur is True:
        return "oui"
    return f"{valeur:g}"


# ---------------- 3. PROGRAMME PRINCIPAL ----------------
if __name__ == "__main__":
    enveloppe, decision = [], []
    for j, (nom, norme, unite, sens) in enumerate(environnements):
        besoin = plus_severe([p[1][j] for p in projets], sens)
        enveloppe.append(besoin)
        if besoin is None:
            decision.append("Pas d'exigence")
        elif couvert(deja_qualifie[1][j], besoin, sens):
            decision.append("Similitude : preuve réutilisée")
        else:
            decision.append("Essai enveloppe à faire")

    sans_mutualisation = sum(v is not None for p in projets for v in p[1])  # un essai par projet et par exigence
    avec_mutualisation = decision.count("Essai enveloppe à faire")

    # note de synthèse
    lignes = ["MATRICE D'ESSAIS : MUTUALISATION DE LA QUALIFICATION (valeurs illustratives)", "=" * 75, "",
              f"Projets : {', '.join(p[0] for p in projets)}",
              f"Preuve existante : {deja_qualifie[0]}", "",
              f"{'Environnement':26s} {'Norme':32s} {'Enveloppe':>10s}  Décision"]
    for j, (nom, norme, unite, sens) in enumerate(environnements):
        lignes.append(f"{nom:26s} {norme:32s} {texte(enveloppe[j]):>6s} {unite:5s} {decision[j]}")
    lignes += ["", f"Essais sans mutualisation (un par projet) : {sans_mutualisation}",
               f"Essais avec mutualisation (enveloppe + similitude) : {avec_mutualisation}",
               f"  -> {sans_mutualisation - avec_mutualisation} essais évités", "",
               "Séquence conseillée sur une seule éprouvette (le feu se fait à part : essai destructif) :"]
    numero = 1
    for j in ordre_essais:
        if decision[j] == "Essai enveloppe à faire":
            lignes.append(f"  {numero}. {environnements[j][0]} ({environnements[j][1]})"
                          f" à {texte(enveloppe[j])} {environnements[j][2]}")
            numero += 1
    lignes += ["", "Critère de réussite de chaque essai : le produit fonctionne avant, pendant et après",
               "l'essai, sans dégradation visible ; ce critère est écrit AVANT l'essai, dans le plan d'essais.",
               "Attention : l'enveloppe est conservative. Si un projet est beaucoup moins sévère que les autres,",
               "le sur-tester peut coûter plus cher (surdimensionnement) qu'un essai dédié."]
    rapport = "\n".join(lignes)
    print(rapport)
    dossier = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats")
    os.makedirs(dossier, exist_ok=True)
    with open(os.path.join(dossier, "matrice_essais.txt"), "w", encoding="utf-8") as fichier:
        fichier.write(rapport)

    # figure : la matrice (une ligne par projet, puis l'enveloppe et la preuve existante)
    noms_lignes = [p[0] for p in projets] + ["ENVELOPPE", deja_qualifie[0]]
    valeurs_lignes = [p[1] for p in projets] + [enveloppe, deja_qualifie[1]]
    couleurs = {"Pas d'exigence": 0, "Similitude : preuve réutilisée": 1, "Essai enveloppe à faire": 2}
    # couleur de chaque colonne selon la décision (gris = rien, vert = similitude, orange = essai)
    fond = np.array([[couleurs[d] for d in decision]] * len(noms_lignes))

    figure, (axe, axe_barres) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [4, 1]})
    axe.imshow(fond, cmap=plt.matplotlib.colors.ListedColormap(["#e0e0e0", "#a8dba8", "#f7c08a"]),
               vmin=0, vmax=2, aspect="auto")
    for i, valeurs in enumerate(valeurs_lignes):
        for j, valeur in enumerate(valeurs):
            gras = "bold" if noms_lignes[i] == "ENVELOPPE" else "normal"
            axe.text(j, i, texte(valeur), ha="center", va="center", fontweight=gras)
    axe.set_xticks(range(len(environnements)),
                   [f"{e[0]}\n[{e[2] or 'oui / non'}]\n{e[1]}" for e in environnements], fontsize=8)
    axe.set_yticks(range(len(noms_lignes)), noms_lignes)
    axe.axhline(len(projets) - 0.5, color="k", lw=2)  # trait entre les projets et l'enveloppe
    axe.set_title("Exigences par projet -> enveloppe -> décision "
                  "(orange : essai à faire, vert : preuve existante réutilisée, gris : pas d'exigence)", fontsize=9)

    axe_barres.bar(["Sans\nmutualisation", "Avec\nmutualisation"], [sans_mutualisation, avec_mutualisation],
                   color=["tab:red", "tab:green"])
    axe_barres.set(title="Nombre d'essais", ylabel="Essais")
    axe_barres.grid(True, axis="y", alpha=0.3)
    figure.suptitle("Matrice d'essais : un référentiel commun pour plusieurs radars (valeurs illustratives)")
    figure.tight_layout()
    figure.savefig(os.path.join(dossier, "matrice_essais.png"), dpi=150)
    plt.show()
