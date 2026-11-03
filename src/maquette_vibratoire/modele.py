"""Maquette : caisson sur silentblocs + carte électronique, sous vibration aléatoire de base.

Modèle à 2 degrés de liberté (translation verticale) :

        x2  ── carte (m2)
              k2, c2       liaison carte / caisson
        x1  ── caisson (m1)
              k1, c1       silentblocs
        y   ── base (porteur), accélération imposée ÿ(t) = DSP d'entrée

Coordonnées relatives à la base z = x - y :
    M z'' + C z' + K z = -M·1·ÿ
En fréquentiel (ÿ = A·e^{iωt}) :
    Z(ω) = (K - ω²M + iωC)⁻¹ (-M·1) · A        → H_rel [m / (m/s²)]
Accélération absolue : Ẍ = A + (-ω²) Z → T(ω) = 1 - ω² H_rel   (transmissibilité, sans unité)

Réponse aléatoire : S_out(f) = |T(f)|² S_in(f),  grandeur RMS = sqrt(∫ S_out df).
L'amortissement est non proportionnel : on calcule les FRF directement (pas de superposition modale).

Unités : SI partout en interne ; la DSP d'entrée est donnée en g²/Hz et convertie avec scipy.constants.g.
"""

import json
import platform
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from numpy.typing import NDArray
from scipy.constants import g as G_SI  # m/s²
from scipy.linalg import eigh

MM_PAR_POUCE = 25.4

type Vecteur = NDArray[np.float64]  # forme (n,)
type Matrice = NDArray[np.float64]  # forme (n, m)
type MatriceComplexe = NDArray[np.complex128]  # forme (n_freq, n_ddl)
type Config = dict[str, Any]  # contenu de config.json
type Resultats = dict[str, Any]  # sortie de analyse_aleatoire


# ---------------------------------------------------------------------------
# Construction du modèle
# ---------------------------------------------------------------------------

def matrices_2ddl(
    m1_kg: float, m2_kg: float, f_susp_Hz: float, zeta_susp: float, f_carte_Hz: float, zeta_carte: float
) -> tuple[Matrice, Matrice, Matrice]:
    """Matrices M, K, C (2x2, SI) du système caisson + carte.

    Convention des paramètres (valeurs « découplées », mesurables séparément) :
      - k1 : raideur des silentblocs telle que l'ensemble (m1+m2) rigide rebondisse à f_susp ;
      - k2 : raideur de la carte telle qu'elle vibre à f_carte sur un caisson rigide ;
      - c_i = 2 ζ_i sqrt(k_i · m_i_associée).
    Les fréquences propres du système couplé diffèrent légèrement de f_susp et f_carte.
    """
    for nom, val in [("m1_kg", m1_kg), ("m2_kg", m2_kg), ("f_susp_Hz", f_susp_Hz), ("f_carte_Hz", f_carte_Hz)]:
        if val <= 0:
            raise ValueError(f"{nom} doit être > 0, reçu {val}")
    for nom, val in [("zeta_susp", zeta_susp), ("zeta_carte", zeta_carte)]:
        if not 0 < val < 1:
            raise ValueError(f"{nom} doit être dans ]0, 1[, reçu {val}")

    m_susp_kg = m1_kg + m2_kg
    k1_N_m = m_susp_kg * (2 * np.pi * f_susp_Hz) ** 2
    k2_N_m = m2_kg * (2 * np.pi * f_carte_Hz) ** 2
    c1_Ns_m = 2 * zeta_susp * np.sqrt(k1_N_m * m_susp_kg)
    c2_Ns_m = 2 * zeta_carte * np.sqrt(k2_N_m * m2_kg)

    M = np.diag([m1_kg, m2_kg])
    K = np.array([[k1_N_m + k2_N_m, -k2_N_m], [-k2_N_m, k2_N_m]])
    C = np.array([[c1_Ns_m + c2_Ns_m, -c2_Ns_m], [-c2_Ns_m, c2_Ns_m]])
    return M, K, C


def modes_propres(M: Matrice, K: Matrice) -> tuple[Vecteur, Matrice]:
    """Fréquences propres [Hz] et déformées normalisées (max |φ| = 1), système non amorti."""
    omega2, phi = eigh(K, M)  # problème généralisé symétrique, valeurs propres croissantes
    assert np.all(omega2 > 0), f"Valeur propre non positive : {omega2}"
    f_Hz = np.sqrt(omega2) / (2 * np.pi)
    phi = phi / np.max(np.abs(phi), axis=0)
    phi = phi * np.sign(phi[0, :])  # signe : composante caisson positive
    return f_Hz, phi


def frf_excitation_base(M: Matrice, K: Matrice, C: Matrice, f_Hz: Vecteur) -> MatriceComplexe:
    """FRF en déplacement relatif par unité d'accélération de base.

    Retourne H_rel de forme (n_freq, n_ddl), en m / (m/s²).
    Le même code sert pour le cas 2 ddl et pour la validation 1 ddl.
    """
    n = M.shape[0]
    assert M.shape == K.shape == C.shape == (n, n)
    assert f_Hz.ndim == 1
    omega = 2 * np.pi * f_Hz
    second_membre = -M @ np.ones(n)  # (n,)
    H = np.empty((f_Hz.size, n), dtype=complex)
    for i, w in enumerate(omega):
        H[i] = np.linalg.solve(K - w**2 * M + 1j * w * C, second_membre)
    return H


def transmissibilite_absolue(H_rel: MatriceComplexe, f_Hz: Vecteur) -> MatriceComplexe:
    """T(ω) = 1 - ω² H_rel : accélération absolue / accélération de base (sans unité)."""
    omega2 = (2 * np.pi * f_Hz) ** 2
    return 1.0 - omega2[:, None] * H_rel


# ---------------------------------------------------------------------------
# DSP d'entrée et réponse aléatoire
# ---------------------------------------------------------------------------

def dsp_trapeze_g2_Hz(
    f_Hz: Vecteur,
    f_min: float,
    f_a: float,
    f_b: float,
    f_max: float,
    niveau_g2_Hz: float,
    pente_montee_dB_oct: float,
    pente_descente_dB_oct: float,
) -> Vecteur:
    """DSP d'accélération en g²/Hz : rampe en dB/octave, palier, rampe ; nulle hors [f_min, f_max]."""
    if not f_min < f_a <= f_b < f_max:
        raise ValueError(f"Fréquences du profil incohérentes : {f_min}, {f_a}, {f_b}, {f_max}")
    # pente en dB/octave → exposant n tel que S ∝ f^n : n = pente / (10 log10 2)
    n_montee = pente_montee_dB_oct / (10 * np.log10(2))
    n_descente = pente_descente_dB_oct / (10 * np.log10(2))
    S = np.zeros_like(f_Hz)
    montee = (f_Hz >= f_min) & (f_Hz < f_a)
    palier = (f_Hz >= f_a) & (f_Hz <= f_b)
    descente = (f_Hz > f_b) & (f_Hz <= f_max)
    S[montee] = niveau_g2_Hz * (f_Hz[montee] / f_a) ** n_montee
    S[palier] = niveau_g2_Hz
    S[descente] = niveau_g2_Hz * (f_Hz[descente] / f_b) ** n_descente
    return S


def rms(dsp: Vecteur, f_Hz: Vecteur) -> float:
    """Valeur RMS d'une DSP unilatérale : sqrt(∫ S df), intégration trapèze."""
    assert dsp.shape == f_Hz.shape
    integrale = np.trapezoid(dsp, f_Hz)
    assert integrale >= 0
    return float(np.sqrt(integrale))


def deplacement_admissible_steinberg_mm(B_mm: float, h_mm: float, L_mm: float, c: float, r: float) -> float:
    """Steinberg : Z_adm = 0.00022 B / (c h r sqrt(L)), toutes longueurs en pouces.

    Si le déplacement 3σ de la carte reste sous Z_adm, Steinberg estime ~2·10⁷ cycles
    avant rupture des brasures (Vibration Analysis for Electronic Equipment, 3e éd.).
    """
    B_in, h_in, L_in = B_mm / MM_PAR_POUCE, h_mm / MM_PAR_POUCE, L_mm / MM_PAR_POUCE
    Z_in = 0.00022 * B_in / (c * h_in * r * np.sqrt(L_in))
    return float(Z_in * MM_PAR_POUCE)


def analyse_aleatoire(M: Matrice, K: Matrice, C: Matrice, f_Hz: Vecteur, S_in_g2_Hz: Vecteur) -> Resultats:
    """Grms caisson/carte et déplacement relatif carte/caisson (RMS) pour une DSP d'entrée."""
    H = frf_excitation_base(M, K, C, f_Hz)
    T = transmissibilite_absolue(H, f_Hz)
    S_acc_g2_Hz = np.abs(T) ** 2 * S_in_g2_Hz[:, None]  # (n_freq, 2)
    S_in_SI = S_in_g2_Hz * G_SI**2  # (m/s²)²/Hz
    d_rel_m_par_ms2 = H[:, 1] - H[:, 0]  # déplacement carte - caisson
    S_drel_m2_Hz = np.abs(d_rel_m_par_ms2) ** 2 * S_in_SI
    return {
        "H": H,
        "T": T,
        "S_acc_g2_Hz": S_acc_g2_Hz,
        "S_drel_m2_Hz": S_drel_m2_Hz,
        "grms_entree": rms(S_in_g2_Hz, f_Hz),
        "grms_caisson": rms(S_acc_g2_Hz[:, 0], f_Hz),
        "grms_carte": rms(S_acc_g2_Hz[:, 1], f_Hz),
        "drel_rms_mm": rms(S_drel_m2_Hz, f_Hz) * 1e3,
    }


# ---------------------------------------------------------------------------
# Validations contre des solutions analytiques
# (appelées à chaque exécution ; elles sont aussi couvertes par les tests pytest)
# ---------------------------------------------------------------------------

def valider_modes_2ddl(M: Matrice, K: Matrice, f_calc_Hz: Vecteur) -> Vecteur:
    """Compare eigh aux racines de det(K - λM) = 0 (formule fermée du 2 ddl)."""
    m1, m2 = M[0, 0], M[1, 1]
    k11, k12, k22 = K[0, 0], K[0, 1], K[1, 1]
    a = m1 * m2
    b = -(m1 * k22 + m2 * k11)
    c = k11 * k22 - k12**2
    lam = np.sort(np.roots([a, b, c]).real)
    f_ana_Hz = np.sqrt(lam) / (2 * np.pi)
    ecart = np.max(np.abs(f_calc_Hz - f_ana_Hz) / f_ana_Hz)
    assert ecart < 1e-9, f"Modes : écart {ecart:.2e} entre eigh et formule fermée"
    return f_ana_Hz


def valider_1ddl_bruit_blanc() -> tuple[float, float, float]:
    """1 ddl sous bruit blanc : ∫|T|² df = π fn (1+4ζ²) / (4ζ) (résultat exact).

    C'est l'équation de Miles corrigée du terme 4ζ² ; elle teste toute la chaîne
    FRF → transmissibilité → intégration avec le même code que le modèle 2 ddl.
    """
    fn_Hz, zeta, m_kg, S0_g2_Hz = 100.0, 0.05, 1.0, 0.01
    k = m_kg * (2 * np.pi * fn_Hz) ** 2
    c = 2 * zeta * np.sqrt(k * m_kg)
    M, K, C = np.array([[m_kg]]), np.array([[k]]), np.array([[c]])
    f_Hz = np.geomspace(0.01, 1e5, 200_000)
    T = transmissibilite_absolue(frf_excitation_base(M, K, C, f_Hz), f_Hz)[:, 0]
    grms_num = rms(np.abs(T) ** 2 * S0_g2_Hz, f_Hz)
    grms_ana = float(np.sqrt(S0_g2_Hz * np.pi * fn_Hz * (1 + 4 * zeta**2) / (4 * zeta)))
    ecart = abs(grms_num - grms_ana) / grms_ana
    assert ecart < 5e-3, f"1 ddl bruit blanc : Grms num {grms_num:.4f} vs analytique {grms_ana:.4f}"
    assert abs(T[0] - 1) < 1e-6, "T(f→0) doit valoir 1 (la masse suit la base)"
    return grms_num, grms_ana, ecart


def valider_convergence_grille(M: Matrice, K: Matrice, C: Matrice, cfg_num: Config, cfg_dsp: Config) -> float:
    """Le Grms carte ne doit pas bouger si on divise le nombre de points par 2."""
    resultats = []
    for n in (cfg_num["n_points"] // 2, cfg_num["n_points"]):
        f = np.geomspace(cfg_num["f_grille_min_Hz"], cfg_num["f_grille_max_Hz"], n)
        resultats.append(analyse_aleatoire(M, K, C, f, dsp_depuis_config(f, cfg_dsp))["grms_carte"])
    ecart = abs(resultats[1] - resultats[0]) / resultats[1]
    assert ecart < 1e-3, f"Grille non convergée : écart {ecart:.2e}"
    return ecart


# ---------------------------------------------------------------------------
# Programme principal
# ---------------------------------------------------------------------------

def dsp_depuis_config(f_Hz: Vecteur, d: Config) -> Vecteur:
    return dsp_trapeze_g2_Hz(
        f_Hz, d["f_min_Hz"], d["f_palier_debut_Hz"], d["f_palier_fin_Hz"], d["f_max_Hz"],
        d["niveau_palier_g2_par_Hz"], d["pente_montante_dB_par_octave"], d["pente_descendante_dB_par_octave"],
    )


def matrices_depuis_config(cfg: Config, zeta_susp: float) -> tuple[Matrice, Matrice, Matrice]:
    return matrices_2ddl(
        cfg["caisson"]["masse_kg"], cfg["carte"]["masse_kg"],
        cfg["caisson"]["frequence_suspension_Hz"], zeta_susp,
        cfg["carte"]["frequence_propre_Hz"], cfg["carte"]["amortissement_ratio"],
    )


def tracer(
    f_Hz: Vecteur,
    S_in: Vecteur,
    res: Resultats,
    f_modes: Vecteur,
    phi: Matrice,
    balayage: list[dict[str, float]],
    dossier: Path,
) -> Path:
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    a = ax[0, 0]
    a.loglog(f_Hz, np.abs(res["T"][:, 0]), label="Caisson")
    a.loglog(f_Hz, np.abs(res["T"][:, 1]), label="Carte")
    a.axhline(1, color="grey", lw=0.8, ls="--")
    for f in f_modes:
        a.axvline(f, color="k", lw=0.6, ls=":")
    a.set(title="Transmissibilité |Ẍ / Ÿ_base|", xlabel="Fréquence [Hz]", ylabel="[-]", xlim=(5, 2000))
    a.legend()
    a.grid(True, which="both", alpha=0.3)

    a = ax[0, 1]
    bande = S_in > 0
    a.loglog(f_Hz[bande], S_in[bande], "k", lw=2, label=f"Entrée ({res['grms_entree']:.2f} Grms)")
    a.loglog(f_Hz[bande], res["S_acc_g2_Hz"][bande, 0], label=f"Caisson ({res['grms_caisson']:.2f} Grms)")
    a.loglog(f_Hz[bande], res["S_acc_g2_Hz"][bande, 1], label=f"Carte ({res['grms_carte']:.2f} Grms)")
    a.set(title="DSP d'accélération", xlabel="Fréquence [Hz]", ylabel="[g²/Hz]")
    a.legend()
    a.grid(True, which="both", alpha=0.3)

    a = ax[1, 0]
    x = np.arange(2)
    for j, f in enumerate(f_modes):
        a.bar(x + (j - 0.5) * 0.35, phi[:, j], width=0.35, label=f"Mode {j + 1} : {f:.1f} Hz")
    a.set_xticks(x, ["Caisson", "Carte"])
    a.axhline(0, color="k", lw=0.8)
    a.set(title="Déformées modales (normalisées)", ylabel="Amplitude [-]")
    a.legend()

    a = ax[1, 1]
    zetas = [b["zeta_susp"] for b in balayage]
    a.plot(zetas, [b["grms_caisson"] for b in balayage], "o-", label="Caisson")
    a.plot(zetas, [b["grms_carte"] for b in balayage], "s-", label="Carte")
    a.axhline(res["grms_entree"], color="k", ls="--", lw=1, label="Entrée")
    a.set(title="Effet de l'amortissement des silentblocs", xlabel="ζ suspension [-]", ylabel="Grms")
    a.legend()
    a.grid(True, alpha=0.3)

    fig.suptitle("Maquette 2 ddl : caisson sur silentblocs + carte électronique (valeurs illustratives)")
    fig.tight_layout()
    chemin = dossier / "synthese.png"
    fig.savefig(chemin, dpi=150)
    plt.close(fig)
    return chemin


def main(chemin_config: str | Path) -> Resultats:
    """Exécute l'étude complète et retourne le dictionnaire écrit dans resultats.json."""
    cfg = json.loads(Path(chemin_config).read_text(encoding="utf-8"))
    dossier = Path(chemin_config).parent / cfg["sortie"]["dossier_figures"]
    dossier.mkdir(exist_ok=True)

    num = cfg["numerique"]
    f_Hz = np.geomspace(num["f_grille_min_Hz"], num["f_grille_max_Hz"], num["n_points"])
    S_in = dsp_depuis_config(f_Hz, cfg["dsp_entree"])

    M, K, C = matrices_depuis_config(cfg, cfg["caisson"]["amortissement_ratio"])
    f_modes, phi = modes_propres(M, K)
    res = analyse_aleatoire(M, K, C, f_Hz, S_in)

    # --- validations (le script s'arrête si l'une échoue)
    valider_modes_2ddl(M, K, f_modes)
    g_num, g_ana, ecart_1ddl = valider_1ddl_bruit_blanc()
    ecart_grille = valider_convergence_grille(M, K, C, num, cfg["dsp_entree"])

    # --- étude paramétrique de l'amortissement des silentblocs
    balayage = []
    for z in cfg["etude_amortissement"]["amortissements_suspension_ratio"]:
        r = analyse_aleatoire(*matrices_depuis_config(cfg, z), f_Hz, S_in)
        balayage.append({"zeta_susp": z, "grms_caisson": r["grms_caisson"], "grms_carte": r["grms_carte"],
                         "drel_3sigma_mm": 3 * r["drel_rms_mm"]})

    st = cfg["critere_steinberg"]
    z_adm_mm = deplacement_admissible_steinberg_mm(
        st["longueur_bord_carte_B_mm"], st["epaisseur_carte_h_mm"], st["longueur_composant_L_mm"],
        st["facteur_composant_c"], st["facteur_position_r"])
    drel_3s_mm = 3 * res["drel_rms_mm"]

    chemin_fig = tracer(f_Hz, S_in, res, f_modes, phi, balayage, dossier)

    # --- rapport console
    T_max = np.abs(res["T"]).max(axis=0)
    print("=== Modes propres (système couplé, non amorti)")
    for j, f in enumerate(f_modes):
        print(f"  Mode {j + 1} : {f:7.2f} Hz   déformée caisson/carte = {phi[0, j]:+.3f} / {phi[1, j]:+.3f}")
    print("=== Transmissibilité maximale")
    print(f"  Caisson : {T_max[0]:.1f}   Carte : {T_max[1]:.1f}")
    print("=== Réponse aléatoire")
    print(f"  Entrée  : {res['grms_entree']:.2f} Grms")
    print(f"  Caisson : {res['grms_caisson']:.2f} Grms")
    print(f"  Carte   : {res['grms_carte']:.2f} Grms   (pic 3σ ≈ {3 * res['grms_carte']:.1f} g)")
    print("=== Critère de Steinberg (carte)")
    print(f"  Déplacement relatif 3σ : {drel_3s_mm:.3f} mm   admissible : {z_adm_mm:.3f} mm   "
          f"marge : {z_adm_mm / drel_3s_mm:.2f}  -> {'OK' if drel_3s_mm <= z_adm_mm else 'NON CONFORME'}")
    print("=== Effet de l'amortissement des silentblocs")
    for b in balayage:
        print(f"  ζ = {b['zeta_susp']:.2f} : caisson {b['grms_caisson']:5.2f} Grms, carte {b['grms_carte']:5.2f} Grms, "
              f"déplacement 3σ {b['drel_3sigma_mm']:.3f} mm")
    print("=== Validations")
    print("  Modes eigh = formule fermée 2 ddl : OK")
    print(f"  1 ddl bruit blanc : Grms num {g_num:.4f} / analytique {g_ana:.4f} (écart {ecart_1ddl:.1e}) : OK")
    print(f"  Convergence grille (n/2 vs n) : écart {ecart_grille:.1e} : OK")
    print(f"Figure : {chemin_fig}")

    resultats = {
        "config": cfg,
        "modes_Hz": f_modes.tolist(),
        "deformees": phi.tolist(),
        "transmissibilite_max": T_max.tolist(),
        "grms": {k: res[k] for k in ("grms_entree", "grms_caisson", "grms_carte")},
        "steinberg": {"drel_3sigma_mm": drel_3s_mm, "admissible_mm": z_adm_mm},
        "balayage_amortissement": balayage,
        "validations": {"ecart_1ddl": ecart_1ddl, "ecart_grille": ecart_grille},
        "environnement": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                          "matplotlib": matplotlib.__version__, "plateforme": platform.platform()},
    }
    (dossier / "resultats.json").write_text(json.dumps(resultats, indent=2, ensure_ascii=False), encoding="utf-8")
    return resultats


def cli() -> None:
    """Point d'entrée : `uv run maquette-vibratoire config.json`."""
    if len(sys.argv) != 2:
        sys.exit("Usage : uv run maquette-vibratoire config.json")
    main(sys.argv[1])


if __name__ == "__main__":
    cli()
