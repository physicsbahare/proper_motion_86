"""Rank alternative JWST/NIRCam epoch pairs for robust PM attempts.

The original selector returns one pair. For very red F444W-selected sources that
can fail when the only same-filter long baseline is in F115W/F150W, this module
provides an auditable ordered list of alternative pairs. Same-filter red bands
remain preferred, but red cross-filter pairs are tried before blue-only pairs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .archive import _group_epochs
from .config import CONFIG, FILTER_WAVELENGTH_NM


def _tier(f1: str, f2: str) -> int:
    w1 = FILTER_WAVELENGTH_NM.get(f1, 0)
    w2 = FILTER_WAVELENGTH_NM.get(f2, 0)
    wmin = min(w1, w2)
    same = f1 == f2
    if same and wmin >= 3560:
        return 0
    if (not same) and wmin >= 3560:
        return 1
    if same and wmin >= 2770:
        return 2
    if (not same) and wmin >= 2770:
        return 3
    if same:
        return 4
    return 5


def ranked_epoch_pairs(inventory: pd.DataFrame, max_pairs: int = 6) -> list[dict]:
    if inventory.empty:
        return []
    collection = inventory.get("obs_collection", pd.Series(index=inventory.index, dtype=object))
    jwst = inventory[collection.astype(str).str.upper().eq("JWST")].copy()
    if "instrument_name" in jwst:
        jwst = jwst[jwst["instrument_name"].astype(str).str.upper().str.contains("NIRCAM")]
    jwst = jwst[jwst["filter_norm"].notna()]
    if jwst.empty:
        return []

    epochs = []
    for filt, grp in jwst.groupby("filter_norm"):
        for e in _group_epochs(grp):
            epochs.append((str(filt), e))

    candidates = []
    for i, (f1, e1) in enumerate(epochs):
        for f2, e2 in epochs[i + 1:]:
            dt = abs(float(e2["mjd"]) - float(e1["mjd"]))
            if dt < CONFIG.min_pm_baseline_days:
                continue
            if e1["mjd"] <= e2["mjd"]:
                ef, early, lf, late = f1, e1, f2, e2
            else:
                ef, early, lf, late = f2, e2, f1, e1
            w1 = FILTER_WAVELENGTH_NM.get(ef, 0)
            w2 = FILTER_WAVELENGTH_NM.get(lf, 0)
            same = ef == lf
            score = (
                _tier(ef, lf),
                -min(w1, w2),
                -(w1 + w2),
                -dt,
            )
            candidates.append((score, {
                "status": "PAIR_FOUND",
                "pair_type": "SAME_FILTER_JWST" if same else "CROSS_FILTER_JWST",
                "filter_early": ef,
                "filter_late": lf,
                "early": early,
                "late": late,
                "baseline_days_inventory": float(dt),
            }))

    candidates.sort(key=lambda x: x[0])
    out = []
    seen = set()
    for _, pair in candidates:
        key = (
            pair["filter_early"], pair["filter_late"],
            round(pair["early"]["mjd"], 4), round(pair["late"]["mjd"], 4),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(pair)
        if len(out) >= max_pairs:
            break
    return out
