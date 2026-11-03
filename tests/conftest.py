import json
import shutil
from pathlib import Path

import pytest

from maquette_vibratoire.modele import Config

RACINE = Path(__file__).parents[1]


@pytest.fixture
def chemin_config() -> Path:
    return RACINE / "config.json"


@pytest.fixture
def config(chemin_config: Path) -> Config:
    return json.loads(chemin_config.read_text(encoding="utf-8"))


@pytest.fixture
def reference() -> dict:
    """Résultats de la version validée, figés avant le passage en projet uv."""
    return json.loads((Path(__file__).parent / "donnees" / "reference.json").read_text(encoding="utf-8"))


@pytest.fixture
def config_temporaire(chemin_config: Path, tmp_path: Path) -> Path:
    """Copie de config.json dans un dossier temporaire : les figures n'écrasent pas celles du projet."""
    copie = tmp_path / "config.json"
    shutil.copy(chemin_config, copie)
    return copie
