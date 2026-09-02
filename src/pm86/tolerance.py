"""Empirical audit of whether a fixed centroid tolerance is method-specific.

The verifier currently discusses a 60 mas position-recovery tolerance. This
module does not assume that 60 mas is correct. It measures two things directly:

1. centroid-method separation: 2-D Gaussian versus center-of-mass centroids;
2. coordinate-frame sensitivity: raw gWCS position versus the same centroid
   after local affine registration.

The independent 282040 centroid measurements can be injected as a reference
calibration so the audit remains informative even before the 45-candidate run
has produced exposure-level astrometry.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_TOLERANCE_MAS = 60.0
SWEEP_MAS = (20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0)


def _num(series):
    return pd.to_numeric(series, errors="coerce")


def _hypot_cols(df: pd.DataFrame, x1: str, y1: str, x2: str, y2: str, scale: float = 1.0):
    needed = [x1, y1, x2, y2]
    if not all(c in df.columns for c in needed):
        return pd.Series(np.nan, index=df.index, dtype=float)
    dx = (_num(df[x1]) - _num(df[x2])) * scale
    dy = (_num(df[y1]) - _num(df[y2])) * scale
    return np.hypot(dx, dy)


def _catalog_radius(df: pd.DataFrame, east: str, north: str, scale: float = 1000.0):
    if east not in df.columns or north not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return np.hypot(_num(df[east]) * scale, _num(df[north]) * scale)


def _small_angle_sep_mas(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    ra1 = pd.to_numeric(ra1_deg, errors="coerce").to_numpy(float)
    dec1 = pd.to_numeric(dec1_deg, errors="coerce").to_numpy(float)
    ra2 = pd.to_numeric(ra2_deg, errors="coerce").to_numpy(float)
    dec2 = pd.to_numeric(dec2_deg, errors="coerce").to_numpy(float)
    mean_dec = np.deg2rad((dec1 + dec2) / 2.0)
    east = (ra1 - ra2) * np.cos(mean_dec) * 3.6e6
    north = (dec1 - dec2) * 3.6e6
    return np.hypot(east, north)


def _load_registered(candidate_dir: Path) -> pd.DataFrame:
    path = candidate_dir / "registered_positions.csv"
    if not path.exists():
        return pd.DataFrame()
    try:
        reg = pd.read_csv(path)
    except Exception:
        return pd.DataFrame()
    if "filename" not in reg.columns:
        return pd.DataFrame()
    keep = [
        c for c in [
            "filename",
            "east_2dg_mas",
            "north_2dg_mas",
            "east_com_mas",
            "north_com_mas",
            "registration_scatter_east_mas",
            "registration_scatter_north_mas",
        ]
        if c in reg.columns
    ]
    return reg[keep].drop_duplicates("filename", keep="first")


def _reference_rows(reference_centroids: str | Path | None) -> pd.DataFrame:
    if reference_centroids is None:
        return pd.DataFrame()
    path = Path(reference_centroids)
    if not path.exists():
        return pd.DataFrame()
    ref = pd.read_csv(path)
    needed = {"ra_gauss", "dec_gauss", "ra_com", "dec_com", "snr_err", "productFilename"}
    if not needed.issubset(ref.columns):
        raise ValueError(f"Reference centroid file is missing columns: {sorted(needed - set(ref.columns))}")

    out = pd.DataFrame(index=ref.index)
    out["source"] = "REFERENCE_282040"
    out["candidate_id"] = pd.to_numeric(ref.get("candidate_id", 282040), errors="coerce").fillna(282040).astype("Int64")
    out["epoch_label"] = "reference"
    out["filter"] = ref.get("filter", "")
    out["filename"] = ref["productFilename"].astype(str)
    out["mjd"] = pd.to_numeric(ref.get("mjd"), errors="coerce")
    out["detector"] = ref.get("detector", "")
    out["clean_snr_err"] = pd.to_numeric(ref["snr_err"], errors="coerce")
    out["raw_snr_err"] = out["clean_snr_err"]
    out["raw_snr_emp"] = np.nan
    out["pixel_scale_mas"] = np.nan
    out["method_sep_2dg_com_mas"] = _small_angle_sep_mas(
        ref["ra_gauss"], ref["dec_gauss"], ref["ra_com"], ref["dec_com"]
    )
    out["registration_shift_2dg_mas"] = np.nan
    out["registration_shift_com_mas"] = np.nan
    out["catalog_offset_2dg_mas"] = np.nan
    out["catalog_offset_com_mas"] = np.nan
    out["registration_scatter_east_mas"] = np.nan
    out["registration_scatter_north_mas"] = np.nan
    return out


def collect_tolerance_rows(
    results_root: str | Path,
    reference_centroids: str | Path | None = None,
) -> pd.DataFrame:
    """Collect one tolerance-audit row per measured exposure plus reference rows."""
    root = Path(results_root)
    rows = []

    for candidate_dir in sorted((root / "candidates").glob("candidate_*")):
        mpath = candidate_dir / "exposure_measurements.csv"
        if not mpath.exists():
            continue
        try:
            m = pd.read_csv(mpath)
        except Exception:
            continue
        if m.empty or "filename" not in m.columns:
            continue

        reg = _load_registered(candidate_dir)
        if len(reg):
            m = m.merge(reg, on="filename", how="left", suffixes=("", "_reg"))

        cid_from_dir = candidate_dir.name.replace("candidate_", "")
        if "candidate_id" not in m.columns:
            m["candidate_id"] = cid_from_dir
        m["source"] = "REMAINING_45"
        m["clean_snr_err"] = _num(m.get("clean_snr_err", pd.Series(np.nan, index=m.index)))
        m["method_sep_2dg_com_mas"] = _hypot_cols(
            m,
            "east_2dg_arcsec_raw_wcs",
            "north_2dg_arcsec_raw_wcs",
            "east_com_arcsec_raw_wcs",
            "north_com_arcsec_raw_wcs",
            scale=1000.0,
        )
        m["catalog_offset_2dg_mas"] = _catalog_radius(
            m, "east_2dg_arcsec_raw_wcs", "north_2dg_arcsec_raw_wcs"
        )
        m["catalog_offset_com_mas"] = _catalog_radius(
            m, "east_com_arcsec_raw_wcs", "north_com_arcsec_raw_wcs"
        )

        if "east_2dg_mas" in m.columns and "north_2dg_mas" in m.columns:
            raw_e = _num(m["east_2dg_arcsec_raw_wcs"]) * 1000.0
            raw_n = _num(m["north_2dg_arcsec_raw_wcs"]) * 1000.0
            m["registration_shift_2dg_mas"] = np.hypot(
                _num(m["east_2dg_mas"]) - raw_e,
                _num(m["north_2dg_mas"]) - raw_n,
            )
        else:
            m["registration_shift_2dg_mas"] = np.nan

        if "east_com_mas" in m.columns and "north_com_mas" in m.columns:
            raw_e = _num(m["east_com_arcsec_raw_wcs"]) * 1000.0
            raw_n = _num(m["north_com_arcsec_raw_wcs"]) * 1000.0
            m["registration_shift_com_mas"] = np.hypot(
                _num(m["east_com_mas"]) - raw_e,
                _num(m["north_com_mas"]) - raw_n,
            )
        else:
            m["registration_shift_com_mas"] = np.nan

        keep = [
            c for c in [
                "source",
                "candidate_id",
                "epoch_label",
                "filter",
                "filename",
                "mjd",
                "detector",
                "clean_snr_err",
                "raw_snr_err",
                "raw_snr_emp",
                "pixel_scale_mas",
                "method_sep_2dg_com_mas",
                "registration_shift_2dg_mas",
                "registration_shift_com_mas",
                "catalog_offset_2dg_mas",
                "catalog_offset_com_mas",
                "registration_scatter_east_mas",
                "registration_scatter_north_mas",
            ]
            if c in m.columns
        ]
        rows.append(m[keep].copy())

    reference = _reference_rows(reference_centroids)
    if len(reference):
        rows.append(reference)

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True, sort=False)
    out["candidate_id"] = pd.to_numeric(out["candidate_id"], errors="coerce").astype("Int64")
    return out


def _pct(values, q):
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, q)) if len(x) else np.nan


def _fraction_le(values, threshold):
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    x = x[np.isfinite(x)]
    return float(np.mean(x <= threshold)) if len(x) else np.nan


def summarize_tolerance(detail: pd.DataFrame, tolerance_mas: float = DEFAULT_TOLERANCE_MAS) -> dict:
    """Summarize the S/N>=3 rows relevant to the verifier's position check."""
    if detail.empty:
        return {
            "tolerance_mas": float(tolerance_mas),
            "eligible_snr3_rows": 0,
            "method_comparable_rows": 0,
            "registration_comparable_rows": 0,
            "assessment": "NOT_EVALUABLE_NO_EXPOSURE_MEASUREMENTS",
        }

    eligible = detail[_num(detail["clean_snr_err"]) >= 3.0].copy()
    method = _num(eligible["method_sep_2dg_com_mas"])
    reg = _num(eligible["registration_shift_2dg_mas"])
    method_finite = method[np.isfinite(method)]
    reg_finite = reg[np.isfinite(reg)]

    method_over = int(np.sum(method_finite > tolerance_mas))
    reg_over = int(np.sum(reg_finite > tolerance_mas))

    if len(method_finite) == 0 and len(reg_finite) == 0:
        assessment = "NOT_EVALUABLE_NO_COMPARABLE_SNR3_ROWS"
    elif method_over == 0 and reg_over == 0:
        assessment = "NO_60MAS_METHOD_OR_FRAME_SENSITIVITY_SEEN"
    else:
        assessment = "60MAS_METHOD_OR_FRAME_SENSITIVITY_DETECTED"

    return {
        "tolerance_mas": float(tolerance_mas),
        "eligible_snr3_rows": int(len(eligible)),
        "eligible_candidates": int(eligible["candidate_id"].nunique()) if len(eligible) else 0,
        "method_comparable_rows": int(len(method_finite)),
        "method_rows_over_tolerance": method_over,
        "method_fraction_within_tolerance": _fraction_le(method_finite, tolerance_mas),
        "method_sep_p50_mas": _pct(method_finite, 50),
        "method_sep_p90_mas": _pct(method_finite, 90),
        "method_sep_p95_mas": _pct(method_finite, 95),
        "method_sep_p99_mas": _pct(method_finite, 99),
        "method_sep_max_mas": _pct(method_finite, 100),
        "registration_comparable_rows": int(len(reg_finite)),
        "registration_rows_over_tolerance": reg_over,
        "registration_fraction_within_tolerance": _fraction_le(reg_finite, tolerance_mas),
        "registration_shift_p50_mas": _pct(reg_finite, 50),
        "registration_shift_p90_mas": _pct(reg_finite, 90),
        "registration_shift_p95_mas": _pct(reg_finite, 95),
        "registration_shift_p99_mas": _pct(reg_finite, 99),
        "registration_shift_max_mas": _pct(reg_finite, 100),
        "assessment": assessment,
    }


def tolerance_sweep(detail: pd.DataFrame, thresholds=SWEEP_MAS) -> pd.DataFrame:
    if detail.empty:
        return pd.DataFrame(
            columns=[
                "tolerance_mas",
                "method_fraction_within",
                "registration_fraction_within",
                "method_comparable_rows",
                "registration_comparable_rows",
            ]
        )
    eligible = detail[_num(detail["clean_snr_err"]) >= 3.0]
    method = _num(eligible["method_sep_2dg_com_mas"])
    method = method[np.isfinite(method)]
    reg = _num(eligible["registration_shift_2dg_mas"])
    reg = reg[np.isfinite(reg)]
    return pd.DataFrame([
        {
            "tolerance_mas": float(t),
            "method_fraction_within": _fraction_le(method, t),
            "registration_fraction_within": _fraction_le(reg, t),
            "method_comparable_rows": int(len(method)),
            "registration_comparable_rows": int(len(reg)),
        }
        for t in thresholds
    ])


def summarize_by_group(detail: pd.DataFrame, group_col: str, tolerance_mas: float = DEFAULT_TOLERANCE_MAS) -> pd.DataFrame:
    if detail.empty or group_col not in detail.columns:
        return pd.DataFrame()
    rows = []
    eligible = detail[_num(detail["clean_snr_err"]) >= 3.0]
    for value, grp in eligible.groupby(group_col, dropna=False):
        s = summarize_tolerance(grp, tolerance_mas)
        s[group_col] = value
        rows.append(s)
    return pd.DataFrame(rows)


def add_tolerance_flags(detail: pd.DataFrame, tolerance_mas: float = DEFAULT_TOLERANCE_MAS) -> pd.DataFrame:
    out = detail.copy()
    out["eligible_snr3"] = _num(out.get("clean_snr_err", pd.Series(np.nan, index=out.index))) >= 3.0
    out["method_sensitive_at_tolerance"] = (
        out["eligible_snr3"]
        & np.isfinite(_num(out.get("method_sep_2dg_com_mas", pd.Series(np.nan, index=out.index))))
        & (_num(out["method_sep_2dg_com_mas"]) > tolerance_mas)
    )
    out["registration_sensitive_at_tolerance"] = (
        out["eligible_snr3"]
        & np.isfinite(_num(out.get("registration_shift_2dg_mas", pd.Series(np.nan, index=out.index))))
        & (_num(out["registration_shift_2dg_mas"]) > tolerance_mas)
    )
    return out
