"""Shared setup for the two figure scripts (house style, paths, small helpers).

The figures are split by narrative, because the two tasks answer different questions and are
read separately:
    scripts/20_figures_status.py -> results/figures/tumour_status/   (S1-S8)
    scripts/21_figures_type.py   -> results/figures/cancer_type/     (T1-T7)
Both read only the tables written by steps 08-19 (and, for panels that plot individual samples,
the local caches). Every figure is SVG + PDF + 400-dpi PNG with `<stem>.csv` of the plotted
values and `<stem>.caption.md`; per-sample values go to `<data_dir>/figure_data/` (D-12).
"""
from __future__ import annotations

import sys

import pandas as pd
from _common import log

from ctsm.config import load_config
from ctsm.figstyle import apply_style

apply_style()
CFG = load_config()
BIO = CFG.results_dir / "biology"
LOCAL = CFG.data_dir / "figure_data"
FIGURES = CFG.results_dir / "figures"
AS_VIEWS = [v for v in CFG.views if v != "gex"]
SCOPE = pd.read_csv(CFG.results_dir / "data" / "task_scope.csv")
CANCERS = list(SCOPE.loc[SCOPE.tumor_status_in_scope, "project"])


def short(region: str) -> str:
    ch, rg = region.split(":")
    a, b = rg.split("-")
    return f"{ch}:{int(a):,}-{int(b):,}"


def run(figs: dict) -> None:
    """Draw the named figures (default: all of this narrative)."""
    for name in sys.argv[1:] or list(figs):
        try:
            figs[name]()
            log(f"{name} written")
        except FileNotFoundError as e:
            log(f"{name} skipped (inputs missing: {e.filename})")
