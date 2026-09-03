"""Local field registration and proper-motion inference."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

from .config import CONFIG
from .measurement import robust_sigma


def _fit_affine(src, dst, min_matches=None, maxiter=8):
    min_matches = min_matches or CONFIG.min_control_matches
    src = np.asarray(src, dtype=float)
    dst = np.asarray(dst, dtype=float)
    good = np.isfinite(src).all(axis=1) & np.isfinite(dst).all(axis=1)
    if good.sum() < min_matches:
        raise RuntimeError(f"only {good.sum()} valid control matches")

    for _ in range(maxiter):
        A = np.column_stack([np.ones(good.sum()), src[good, 0], src[good, 1]])
        coef = np.linalg.lstsq(A, dst[good], rcond=None)[0]
        Aall = np.column_stack([np.ones(len(src)), src[:, 0], src[:, 1]])
        pred = Aall @ coef
        resid = dst - pred
        rr = np.hypot(resid[:, 0], resid[:, 1])
        med = np.nanmedian(rr[good])
        sig = robust_sigma(rr[good])
        if not np.isfinite(sig) or sig <= 0:
            break
        threshold = max(med + 4.0 * sig, 0.003)
        new_good = good & (rr < threshold)
        if new_good.sum() < min_matches or np.array_equal(new_good, good):
            break
        good = new_good

    A = np.column_stack([np.ones(good.sum()), src[good, 0], src[good, 1]])
    coef = np.linalg.lstsq(A, dst[good], rcond=None)[0]
    Aall = np.column_stack([np.ones(len(src)), src[:, 0], src[:, 1]])
    resid = dst - Aall @ coef
    sx = robust_sigma(resid[good, 0])
    sy = robust_sigma(resid[good, 1])
    return coef, good, resid, float(sx), float(sy)


def _apply_affine(coef, east, north):
    return np.array([1.0, float(east), float(north)]) @ np.asarray(coef)


def _match_controls(current, reference):
    if len(current) < CONFIG.min_control_matches or len(reference) < CONFIG.min_control_matches:
        raise RuntimeError("too few controls before matching")
    src = current[["east_arcsec_raw_wcs", "north_arcsec_raw_wcs"]].to_numpy(float)
    ref = reference[["east_arcsec_raw_wcs", "north_arcsec_raw_wcs"]].to_numpy(float)
    tree = cKDTree(ref)
    dist, idx = tree.query(src, k=1)
    m = pd.DataFrame({"cidx": np.arange(len(src)), "ridx": idx, "sep": dist})
    m = (
        m[m["sep"] < CONFIG.match_radius_arcsec]
        .sort_values("sep")
        .drop_duplicates("ridx", keep="first")
        .reset_index(drop=True)
    )
    if len(m) < CONFIG.min_control_matches:
        raise RuntimeError(f"only {len(m)} unique matched controls")
    return m["cidx"].to_numpy(), m["ridx"].to_numpy()


def register_exposures(measurements: pd.DataFrame, controls: pd.DataFrame):
    """Register every usable exposure into an early-epoch reference frame.

    ``astrometric_snr`` is the authoritative first26 target-significance column
    when present.  It can represent a quality-controlled forced Gaussian fit for
    a faint source whose aperture ``clean_snr_err`` is below the nominal gate.
    Other workflows that do not provide this column retain the original aperture
    S/N behaviour.
    """
    if measurements.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    if "astrometric_snr" in measurements.columns:
        target_snr = pd.to_numeric(measurements["astrometric_snr"], errors="coerce")
    else:
        target_snr = pd.to_numeric(measurements["clean_snr_err"], errors="coerce")

    usable = measurements[
        np.isfinite(measurements["east_2dg_arcsec_raw_wcs"])
        & np.isfinite(measurements["north_2dg_arcsec_raw_wcs"])
        & (target_snr >= CONFIG.target_min_snr_astrometry)
    ].copy()
    if usable.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    counts = controls.groupby("filename").size() if len(controls) else pd.Series(dtype=int)
    usable["n_controls"] = usable["filename"].map(counts).fillna(0)
    early = usable[usable["epoch_label"].eq("early")]
    pool = early if len(early) else usable
    reference_row = pool.sort_values(["n_controls", "mjd"], ascending=[False, True]).iloc[0]
    reference_file = str(reference_row["filename"])
    reference_controls = controls[controls["filename"].eq(reference_file)].reset_index(drop=True)

    registered = []
    quality = []
    matched_rows = []

    for _, row in usable.sort_values("mjd").iterrows():
        filename = str(row["filename"])
        current_controls = controls[controls["filename"].eq(filename)].reset_index(drop=True)

        if filename == reference_file:
            coef = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
            n_match = len(reference_controls)
            n_used = n_match
            sx = sy = 0.0
            method = "REFERENCE_IDENTITY"
            color_r_e = color_r_n = color_p_e = color_p_n = np.nan
        else:
            cidx, ridx = _match_controls(current_controls, reference_controls)
            cur = current_controls.iloc[cidx]
            ref = reference_controls.iloc[ridx]
            src = cur[["east_arcsec_raw_wcs", "north_arcsec_raw_wcs"]].to_numpy(float)
            dst = ref[["east_arcsec_raw_wcs", "north_arcsec_raw_wcs"]].to_numpy(float)
            coef, good, resid, sx, sy = _fit_affine(src, dst)
            n_match = len(src)
            n_used = int(good.sum())
            method = "LOCAL_AFFINE"

            color_r_e = color_r_n = color_p_e = color_p_n = np.nan
            if str(row["filter"]) != str(reference_row["filter"]):
                fcur = cur["flux"].to_numpy(float)[good]
                fref = ref["flux"].to_numpy(float)[good]
                ok = (fcur > 0) & (fref > 0) & np.isfinite(fcur) & np.isfinite(fref)
                if ok.sum() >= 6:
                    color = -2.5 * np.log10(fcur[ok] / fref[ok])
                    try:
                        color_r_e, color_p_e = pearsonr(color, resid[good, 0][ok])
                        color_r_n, color_p_n = pearsonr(color, resid[good, 1][ok])
                    except Exception:
                        pass

            for k, (ci, ri) in enumerate(zip(cidx, ridx)):
                matched_rows.append({
                    "filename": filename,
                    "reference_file": reference_file,
                    "current_control_index": int(ci),
                    "reference_control_index": int(ri),
                    "used": bool(good[k]),
                    "resid_east_arcsec": float(resid[k, 0]),
                    "resid_north_arcsec": float(resid[k, 1]),
                })

        e2, n2 = _apply_affine(coef, row["east_2dg_arcsec_raw_wcs"], row["north_2dg_arcsec_raw_wcs"])
        ec, nc = (np.nan, np.nan)
        if np.isfinite(row["east_com_arcsec_raw_wcs"]) and np.isfinite(row["north_com_arcsec_raw_wcs"]):
            ec, nc = _apply_affine(coef, row["east_com_arcsec_raw_wcs"], row["north_com_arcsec_raw_wcs"])

        registered.append({
            **row.to_dict(),
            "reference_file": reference_file,
            "east_2dg_mas": float(e2 * 1000.0),
            "north_2dg_mas": float(n2 * 1000.0),
            "east_com_mas": float(ec * 1000.0) if np.isfinite(ec) else np.nan,
            "north_com_mas": float(nc * 1000.0) if np.isfinite(nc) else np.nan,
            "registration_scatter_east_mas": float(sx * 1000.0),
            "registration_scatter_north_mas": float(sy * 1000.0),
        })
        quality.append({
            "filename": filename,
            "filter": row["filter"],
            "mjd": row["mjd"],
            "reference_file": reference_file,
            "method": method,
            "n_match": n_match,
            "n_used": n_used,
            "scatter_east_mas": sx * 1000.0,
            "scatter_north_mas": sy * 1000.0,
            "color_corr_east_r": color_r_e,
            "color_corr_east_p": color_p_e,
            "color_corr_north_r": color_r_n,
            "color_corr_north_p": color_p_n,
        })

    return pd.DataFrame(registered), pd.DataFrame(quality), pd.DataFrame(matched_rows)


def _weighted_epoch(values, sigmas):
    v = np.asarray(values, dtype=float)
    s = np.asarray(sigmas, dtype=float)
    ok = np.isfinite(v) & np.isfinite(s) & (s > 0)
    if ok.sum() == 0:
        return np.nan, np.nan
    w = 1.0 / s[ok] ** 2
    mean = np.sum(w * v[ok]) / np.sum(w)
    formal = 1.0 / np.sqrt(np.sum(w))
    empirical = np.std(v[ok], ddof=1) / np.sqrt(ok.sum()) if ok.sum() > 1 else 0.0
    return float(mean), float(max(formal, empirical))


def _position_errors(reg: pd.DataFrame):
    # Use the local gWCS pixel scale measured at the target in each exposure.
    # Nominal NIRCam LW/SW values are only a fallback if gWCS evaluation fails.
    nominal = np.where(
        reg["filter"].isin(["F277W", "F356W", "F410M", "F444W"]),
        63.0,
        31.0,
    )
    if "pixel_scale_mas" in reg.columns:
        measured = pd.to_numeric(reg["pixel_scale_mas"], errors="coerce").to_numpy(float)
        scale = np.where(np.isfinite(measured) & (measured > 0), measured, nominal)
    else:
        scale = nominal

    sx = pd.to_numeric(reg["sigma_x_pix"], errors="coerce").to_numpy(float)
    sy = pd.to_numeric(reg["sigma_y_pix"], errors="coerce").to_numpy(float)
    sx = np.where(np.isfinite(sx) & (sx > 0), sx * scale, 8.0)
    sy = np.where(np.isfinite(sy) & (sy > 0), sy * scale, 8.0)

    rx = pd.to_numeric(reg["registration_scatter_east_mas"], errors="coerce").to_numpy(float)
    ry = pd.to_numeric(reg["registration_scatter_north_mas"], errors="coerce").to_numpy(float)
    positive_rx = rx[np.isfinite(rx) & (rx > 0)]
    positive_ry = ry[np.isfinite(ry) & (ry > 0)]
    floor_x = float(np.nanmedian(positive_rx)) if len(positive_rx) else 5.0
    floor_y = float(np.nanmedian(positive_ry)) if len(positive_ry) else 5.0
    rx = np.where(np.isfinite(rx) & (rx > 0), rx, floor_x)
    ry = np.where(np.isfinite(ry) & (ry > 0), ry, floor_y)
    return np.sqrt(sx**2 + rx**2), np.sqrt(sy**2 + ry**2)


def _fit_method(reg, east_col, north_col, sigma_e, sigma_n):
    early = reg[reg["epoch_label"].eq("early")].copy()
    late = reg[reg["epoch_label"].eq("late")].copy()
    if len(early) == 0 or len(late) == 0:
        return None, []

    early = early.assign(_se=sigma_e[early.index], _sn=sigma_n[early.index])
    late = late.assign(_se=sigma_e[late.index], _sn=sigma_n[late.index])

    ee, see = _weighted_epoch(early[east_col], early["_se"])
    en, sen = _weighted_epoch(early[north_col], early["_sn"])
    le, sle = _weighted_epoch(late[east_col], late["_se"])
    ln, sln = _weighted_epoch(late[north_col], late["_sn"])
    if not all(np.isfinite(x) for x in [ee, see, en, sen, le, sle, ln, sln]):
        return None, []

    baseline_days = float(late["mjd"].mean() - early["mjd"].mean())
    if baseline_days < CONFIG.min_pm_baseline_days:
        return None, []
    dt = baseline_days / 365.25
    de, dn = le - ee, ln - en
    sde, sdn = np.hypot(see, sle), np.hypot(sen, sln)
    result = {
        "baseline_days": baseline_days,
        "baseline_year": dt,
        "delta_east_mas": de,
        "delta_north_mas": dn,
        "sigma_delta_east_mas": sde,
        "sigma_delta_north_mas": sdn,
        "mu_alpha_cosdec_masyr": de / dt,
        "mu_delta_masyr": dn / dt,
        "sigma_mu_alpha_cosdec_masyr": sde / dt,
        "sigma_mu_delta_masyr": sdn / dt,
        "significance_2d": float(np.hypot(de / sde, dn / sdn)),
    }

    pairs = []
    for _, a in early.iterrows():
        for _, b in late.iterrows():
            dty = (b["mjd"] - a["mjd"]) / 365.25
            if dty <= 0:
                continue
            pde = b[east_col] - a[east_col]
            pdn = b[north_col] - a[north_col]
            pairs.append({
                "early_file": a["filename"],
                "late_file": b["filename"],
                "baseline_year": dty,
                "delta_east_mas": pde,
                "delta_north_mas": pdn,
                "mu_alpha_cosdec_masyr": pde / dty,
                "mu_delta_masyr": pdn / dty,
            })
    return result, pairs


def infer_proper_motion(registered: pd.DataFrame, quality: pd.DataFrame, pair_type: str):
    if registered.empty:
        return {"pm_status": "INSUFFICIENT_DATA", "reason": "no registered target positions"}, pd.DataFrame()

    reg = registered.reset_index(drop=True).copy()
    sigma_e, sigma_n = _position_errors(reg)
    result_2dg, pairs_2dg = _fit_method(reg, "east_2dg_mas", "north_2dg_mas", sigma_e, sigma_n)
    result_com, pairs_com = _fit_method(reg, "east_com_mas", "north_com_mas", sigma_e, sigma_n)

    if result_2dg is None:
        return {"pm_status": "INSUFFICIENT_DATA", "reason": "two usable registered epochs not available"}, pd.DataFrame()

    pairwise = pd.DataFrame([{**p, "method": "2DG"} for p in pairs_2dg] + [{**p, "method": "COM"} for p in pairs_com])

    p2 = pairwise[pairwise["method"].eq("2DG")]
    if len(p2):
        mean_e = result_2dg["mu_alpha_cosdec_masyr"]
        mean_n = result_2dg["mu_delta_masyr"]
        norm = np.hypot(mean_e, mean_n)
        if norm > 0:
            dot = (p2["mu_alpha_cosdec_masyr"] * mean_e + p2["mu_delta_masyr"] * mean_n) / norm
            direction_fraction = float(np.mean(dot > 0))
        else:
            direction_fraction = np.nan
    else:
        direction_fraction = np.nan

    method_disagreement = False
    if result_com is not None:
        dmu = np.hypot(
            result_2dg["mu_alpha_cosdec_masyr"] - result_com["mu_alpha_cosdec_masyr"],
            result_2dg["mu_delta_masyr"] - result_com["mu_delta_masyr"],
        )
        comb = np.hypot(
            np.hypot(result_2dg["sigma_mu_alpha_cosdec_masyr"], result_com["sigma_mu_alpha_cosdec_masyr"]),
            np.hypot(result_2dg["sigma_mu_delta_masyr"], result_com["sigma_mu_delta_masyr"]),
        )
        method_disagreement = bool(dmu > max(10.0, 3.0 * comb))
    else:
        dmu = np.nan

    color_systematic = False
    if pair_type == "CROSS_FILTER_JWST" and len(quality):
        for axis in ["east", "north"]:
            r = pd.to_numeric(quality[f"color_corr_{axis}_r"], errors="coerce")
            p = pd.to_numeric(quality[f"color_corr_{axis}_p"], errors="coerce")
            if ((r.abs() >= CONFIG.cross_filter_color_corr_limit) & (p < 0.05)).any():
                color_systematic = True

    sig = result_2dg["significance_2d"]
    max_sigma = max(
        result_2dg["sigma_mu_alpha_cosdec_masyr"],
        result_2dg["sigma_mu_delta_masyr"],
    )

    if (
        sig >= CONFIG.moving_significance
        and (not np.isfinite(direction_fraction) or direction_fraction >= CONFIG.min_pair_direction_fraction)
        and not method_disagreement
        and not color_systematic
    ):
        classification = "MOVING"
    elif sig < CONFIG.zero_motion_significance and max_sigma <= CONFIG.max_stationary_sigma_masyr:
        classification = "CONSISTENT_WITH_ZERO"
    elif method_disagreement or color_systematic:
        classification = "AMBIGUOUS_SYSTEMATICS"
    else:
        classification = "AMBIGUOUS"

    evidence_grade = "A_SAME_FILTER_JWST" if pair_type == "SAME_FILTER_JWST" else "B_CROSS_FILTER_JWST"
    summary = {
        "pm_status": "PM_MEASURED",
        "classification": classification,
        "evidence_grade": evidence_grade,
        "pair_type": pair_type,
        "pair_direction_fraction": direction_fraction,
        "centroid_method_disagreement": method_disagreement,
        "centroid_method_delta_masyr": dmu,
        "cross_filter_color_systematic_flag": color_systematic,
        **result_2dg,
    }
    if result_com is not None:
        summary.update({f"com_{k}": v for k, v in result_com.items()})
    return summary, pairwise
