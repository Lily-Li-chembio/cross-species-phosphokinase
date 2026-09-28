"""Shared paths and helpers for the numbered pipeline scripts."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))  # run the scripts without installing the package

DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "docs" / "figures"

# Default input locations (see data/README.md for where to download each file)
KINOME_MATRICES = DATA / "kinome_atlas" / "ST2.xlsx"      # Johnson et al. 2023, Supp. Table 2
KINOME_REFERENCE = DATA / "kinome_atlas" / "ST3.csv"      # Johnson et al. 2023, Supp. Table 3
MOUSE_ATLAS = DATA / "mouse_atlas" / "ST2.xlsx"           # Giansanti et al. 2022, Supp. Table 2
MOUSE_UNIPROT = DATA / "uniprot" / "mouse.tsv"
HUMAN_UNIPROT = DATA / "uniprot" / "human.tsv"


def out(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
