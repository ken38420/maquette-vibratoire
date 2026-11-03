"""Tests de la maquette : solutions analytiques, cas limites, erreurs attendues, non-régression.

Chaque test compare le code à un résultat connu par ailleurs (formule exacte ou calcul à la main),
jamais à lui-même, sauf le test de non-régression qui le compare à la version validée.
"""

from pathlib import Path

import numpy as np
import pytest

from maquette_vibratoire.modele import (
    Config,
    cli,
    deplacement_admissible_steinberg_mm,
    dsp_depuis_config,
    dsp_trapeze_g2_Hz,
    frf_excitation_base,
    main,
    matrices_2ddl,
    matrices_depuis_config,
    modes_propres,
    rms,
    transmissibilite_absolue,
    valider_1ddl_bruit_blanc,
    valider_convergence_grille,
    valider_modes_2ddl,
)


def transmissibilite_1ddl(fn_Hz: float, zeta: float, f_Hz: np.ndarray) -> np.ndarray:
    """|T| d'un 1 ddl calculé par le code du modèle (masse 1 kg, sans perte de généralité)."""
    k = (2 * np.pi * fn_Hz) ** 2
    c = 2 * zeta * np.sqrt(k)
    M, K, C = np.array([[1.0]]), np.array([[k]]), np.array([[c]])
    return np.abs(transmissibilite_absolue(frf_excitation_base(M, K, C, f_Hz), f_Hz)[:, 0])


# ---------------------------------------------------------------------------
# Paramètres invalides : le code doit s'arrêter, pas produire un résultat faux
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "parametres, message",
    [
        ((0.0, 0.3, 25, 0.1, 200, 0.03), "m1_kg"),
        ((25, -0.3, 25, 0.1, 200, 0.03), "m2_kg"),
        ((25, 0.3, 0.0, 0.1, 200, 0.03), "f_susp_Hz"),
        ((25, 0.3, 25, 0.0, 200, 0.03), "zeta_susp"),
        ((25, 0.3, 25, 0.1, 200, 1.0), "zeta_carte"),
    ],
)
def test_matrices_2ddl_refuse_parametres_invalides(parametres: tuple, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        matrices_2ddl(*parametres)


def test_dsp_refuse_profil_incoherent() -> None:
    with pytest.raises(ValueError, match="incohérentes"):
        dsp_trapeze_g2_Hz(np.array([100.0]), 20, 400, 300, 2000, 0.04, 3, -3)


def test_cli_sans_argument_affiche_usage(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.argv", ["maquette-vibratoire"])
    with pytest.raises(SystemExit, match="Usage"):
        cli()


# ---------------------------------------------------------------------------
# Solutions analytiques et cas limites
# ---------------------------------------------------------------------------

def test_modes_egaux_a_la_formule_fermee(config: Config) -> None:
    M, K, _ = matrices_depuis_config(config, config["caisson"]["amortissement_ratio"])
    f_Hz, _ = modes_propres(M, K)
    f_ana_Hz = valider_modes_2ddl(M, K, f_Hz)
    np.testing.assert_allclose(f_Hz, f_ana_Hz, rtol=1e-9)


def test_1ddl_bruit_blanc_egal_a_la_solution_exacte() -> None:
    grms_num, grms_ana, _ = valider_1ddl_bruit_blanc()
    assert grms_num == pytest.approx(grms_ana, rel=5e-3)


@pytest.mark.parametrize("zeta", [0.02, 0.05, 0.1, 0.3])
def test_1ddl_transmissibilite_a_la_resonance(zeta: float) -> None:
    """À f = fn : |T| = sqrt(1 + 4ζ²) / (2ζ) (formule exacte)."""
    fn_Hz = 100.0
    T = transmissibilite_1ddl(fn_Hz, zeta, np.array([fn_Hz]))
    assert T[0] == pytest.approx(np.sqrt(1 + 4 * zeta**2) / (2 * zeta), rel=1e-12)


@pytest.mark.parametrize("zeta", [0.02, 0.05, 0.1, 0.3])
def test_1ddl_isolation_commence_a_racine_de_deux(zeta: float) -> None:
    """À f = √2·fn : |T| = 1 quel que soit l'amortissement (début de l'isolation)."""
    fn_Hz = 100.0
    T = transmissibilite_1ddl(fn_Hz, zeta, np.array([np.sqrt(2) * fn_Hz]))
    assert T[0] == pytest.approx(1.0, rel=1e-12)


def test_2ddl_basse_frequence_suit_la_base(config: Config) -> None:
    """f → 0 : caisson et carte ont la même accélération que la base (T = 1)."""
    M, K, C = matrices_depuis_config(config, config["caisson"]["amortissement_ratio"])
    f_Hz = np.array([1e-3])
    T = transmissibilite_absolue(frf_excitation_base(M, K, C, f_Hz), f_Hz)
    np.testing.assert_allclose(np.abs(T), 1.0, atol=1e-8)


def test_dsp_pente_3dB_par_octave_et_palier() -> None:
    """+3 dB/octave : une octave sous le palier, la DSP vaut niveau × 10^(-0.3). Nulle hors bande."""
    niveau = 0.04
    f_Hz = np.array([10.0, 40.0, 100.0, 700.0, 3000.0])
    S = dsp_trapeze_g2_Hz(f_Hz, 20, 80, 350, 2000, niveau, 3, -3)
    attendu = [0.0, niveau * 10 ** (-0.3), niveau, niveau * 10 ** (-0.3), 0.0]
    np.testing.assert_allclose(S, attendu, rtol=1e-12)


def test_rms_dsp_constante() -> None:
    """DSP constante S0 sur [a, b] : RMS = sqrt(S0 (b - a))."""
    f_Hz = np.linspace(20.0, 2000.0, 101)
    assert rms(np.full_like(f_Hz, 0.04), f_Hz) == pytest.approx(np.sqrt(0.04 * 1980.0), rel=1e-12)


def test_steinberg_toutes_dimensions_un_pouce() -> None:
    """B = h = L = 1 pouce, c = r = 1 : Z_adm = 0.00022 pouce = 0.005588 mm."""
    assert deplacement_admissible_steinberg_mm(25.4, 25.4, 25.4, 1.0, 1.0) == pytest.approx(0.00022 * 25.4)


def test_grille_frequentielle_convergee(config: Config) -> None:
    M, K, C = matrices_depuis_config(config, config["caisson"]["amortissement_ratio"])
    assert valider_convergence_grille(M, K, C, config["numerique"], config["dsp_entree"]) < 1e-3


def test_dsp_entree_grms_coherent_avec_integrale_manuelle(config: Config) -> None:
    """Grms d'entrée recalculé par intégration analytique des trois segments en loi puissance."""
    d = config["dsp_entree"]
    f_min, f_a, f_b, f_max = d["f_min_Hz"], d["f_palier_debut_Hz"], d["f_palier_fin_Hz"], d["f_max_Hz"]
    S0 = d["niveau_palier_g2_par_Hz"]
    n_m = d["pente_montante_dB_par_octave"] / (10 * np.log10(2))
    n_d = d["pente_descendante_dB_par_octave"] / (10 * np.log10(2))

    def integrale_loi_puissance(f_ref: float, n: float, f1: float, f2: float) -> float:
        return S0 * f_ref / (n + 1) * ((f2 / f_ref) ** (n + 1) - (f1 / f_ref) ** (n + 1))

    aire = (integrale_loi_puissance(f_a, n_m, f_min, f_a) + S0 * (f_b - f_a)
            + integrale_loi_puissance(f_b, n_d, f_b, f_max))
    f_Hz = np.geomspace(1.0, 3000.0, 40_000)
    assert rms(dsp_depuis_config(f_Hz, d), f_Hz) == pytest.approx(np.sqrt(aire), rel=1e-4)


# ---------------------------------------------------------------------------
# Non-régression : l'étude complète redonne les résultats de la version validée
# ---------------------------------------------------------------------------

def test_non_regression_etude_complete(config_temporaire: Path, reference: dict) -> None:
    res = main(config_temporaire)

    np.testing.assert_allclose(res["modes_Hz"], reference["modes_Hz"], rtol=1e-9)
    np.testing.assert_allclose(res["transmissibilite_max"], reference["transmissibilite_max"], rtol=1e-9)
    for cle, valeur in reference["grms"].items():
        assert res["grms"][cle] == pytest.approx(valeur, rel=1e-9), cle
    for cle, valeur in reference["steinberg"].items():
        assert res["steinberg"][cle] == pytest.approx(valeur, rel=1e-9), cle
    assert len(res["balayage_amortissement"]) == len(reference["balayage_amortissement"])
    for obtenu, attendu in zip(res["balayage_amortissement"], reference["balayage_amortissement"], strict=True):
        for cle, valeur in attendu.items():
            assert obtenu[cle] == pytest.approx(valeur, rel=1e-9), cle

    assert (config_temporaire.parent / "figures" / "synthese.png").stat().st_size > 0
