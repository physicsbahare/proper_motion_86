#!/usr/bin/env python3
"""Galaxy-only leave-one-out cross-filter registration diagnostic.

Read the user's LOCAL COSMOS2025-matched controls CSV, not a GitHub
copy of the large catalog. Use only photometrically-classified galaxy
references. Keep the original manuscript measurements unchanged.
This is a representative-exposure-pair diagnostic, not full-epoch
systematic covariance, stationary verification, or a formal PM significance.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


REQUIRED = {
    "candidate_id", "galaxy_photometric_candidate", "catalog_position_match",
    "early_east_arcsec", "early_north_arcsec",
    "late_east_arcsec", "late_north_arcsec",
    "snr_early", "snr_late", "retained_for_covariance",
    "east_mas", "north_mas"
}


def positive_bool(values):
    if values.dtype == bool:
        return values
    return values.astype(str).str.lower().isin(["true", "1", "yes"])


def robust_registration(src: np.ndarray, dst: np.ndarray, min_fit: int = 6):
    good = np.isfinite(src).all(axis=1) & np.isfinite(dst).all(axis=1)
    if good.sum() < min_fit:
        raise ValueError("INSUFFICIENT_FIT_CONTROLS")
    for _ in range(8):
        a = np.c_[np.ones(int(good.sum())), src[good]]
        if np.linalg.matrix_rank(a) < 3:
            raise ValueError("DEGENERATE_FIT_GEOMETRY")
        coef = np.linalg.lstsq(a, dst[good], rcond=None)[0]
        err = dst - np.c_[np.ones(len(src)), src] @ coef
        radial = np.linalg.norm(err, axis=1)
        central = float(np.median(radial[good]))
        mad = float(1.4826 * np.median(np.abs(radial[good] - central)))
        threshold = max(central + 4 * mad, .003)
        new_good = good & (radial <= threshold)
        if new_good.sum() < min_fit or np.array_equal(new_good, good):
            break
        good = new_good
    a = np.c_[np.ones(int(good.sum())), src[good]]
    if np.linalg.matrix_rank(a) < 3:
        raise ValueError("DEGENERATE_FIT_GEOMETRY")
    coef = np.linalg.lstsq(a, dst[good], rcond=None)[0]
    return coef, int(good.sum())


def galaxy_loo(group: pd.DataFrame, min_fit: int = 6):
    selected = group.loc[
        positive_bool(group["galaxy_photometric_candidate"]) &
        positive_bool(group["catalog_position_match"]) &
        (pd.to_numeric(group["snr_early"], errors="coerce") >= 10) &
        (pd.to_numeric(group["snr_late"], errors="coerce") >= 10)
    ].copy()
    n = len(selected)
    result = {
        "n_all_controls": int(len(group)), "n_photometric_galaxies": n,
        "n_original_covariance_retained": int(positive_bool(selected["retained_for_covariance"]).sum()),
        "spectroscopically_confirmed_stationary": 0,
        "full_epoch_systematics_calibrated": False,
        "selection": "catalog_match_and_LePHARE_type0_and_SNR>=10",
    }
    if n < min_fit + 1:
        return pd.DataFrame(), {**result,"status":"INSUFFICIENT_GALAXY_CONTROLS_FOR_LOO","min_required":min_fit+1}
    src = selected[["late_east_arcsec","late_north_arcsec"]].to_numpy(float)
    dst = selected[["early_east_arcsec","early_north_arcsec"]].to_numpy(float)
    rows=[]
    for i in range(n):
        keep = np.arange(n) != i
        try:
            coef,fit_n=robust_registration(src[keep],dst[keep],min_fit)
            xy=np.r_[1,src[i]]@coef
            resid=(xy-dst[i])*1000
            rows.append({
                "candidate_id":int(selected.iloc[i]["candidate_id"]),
                "local_control_index":selected.iloc[i].get("local_control_index",i),
                "catalog_id":selected.iloc[i].get("id"),
                "catalog_match_sep_arcsec":selected.iloc[i].get("catalog_match_sep_arcsec"),
                "lephare_zphot":selected.iloc[i].get("lephare_zphot"),
                "galaxy_type":selected.iloc[i].get("lephare_type"),
                "original_retained_for_covariance":bool(
                    positive_bool(pd.Series([selected.iloc[i]["retained_for_covariance"]])).iloc[0]),
                "n_fit_galaxies":fit_n,
                "galaxy_only_east_mas":float(resid[0]),
                "galaxy_only_north_mas":float(resid[1]),
                "original_mixed_east_mas":float(selected.iloc[i]["east_mas"]),
                "original_mixed_north_mas":float(selected.iloc[i]["north_mas"]),
            })
        except (ValueError, np.linalg.LinAlgError):
            pass
    out=pd.DataFrame(rows)
    result["n_loo_success"]=int(len(out))
    if len(out)<min_fit+1:
        return out,{**result,"status":"INSUFFICIENT_SUCCESSFUL_GALAXY_LOO"}
    xy=out[["galaxy_only_east_mas","galaxy_only_north_mas"]].to_numpy(float)
    mean=np.median(xy,axis=0)
    cov=np.cov(xy,rowvar=False,ddof=1)
    result.update({
        "status":"PRELIMINARY_GALAXY_ONLY_PAIR_NULL",
        "median_east_north_mas":mean.tolist(),
        "scatter_east_north_mas":np.sqrt(np.diag(cov)).tolist(),
        "covariance_mas2":cov.tolist(),
        "n_loo_success":len(out),
        "warning":"LePHARE galaxies are photometric candidates, not independently confirmed stationary sources; results mix centroid noise and registration. Do not add entire covariance to frozen errors or call output Gaussian sigma."
    })
    return out,result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",default="results_stationary_galaxy_controls/cosmos2025_matched_field_controls.csv")
    ap.add_argument("--out",default="results_stationary_galaxy_controls")
    args=ap.parse_args()
    data=pd.read_csv(args.input)
    missing=REQUIRED-set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}. Re-run updated COSMOS2025 crossmatch before this step.")
    if not data["lephare_type"].notna().any():
        raise ValueError("Missing LePHARE classifications. Re-run COSMOS2025 match on latest review branch.")
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    records=[];summaries={}
    for cid,group in data.groupby("candidate_id"):
        residual,summary=galaxy_loo(group)
        summaries[str(cid)]=summary
        if not residual.empty:
            records.append(residual)
    if records:
        pd.concat(records,ignore_index=True).to_csv(out/"galaxy_only_loo_residuals.csv",index=False)
    result={
        "input":str(args.input), "method":"galaxy-only affine leave-one-out using original representative detector exposure pair",
        "candidate_results":summaries,
        "publication_status":"DIAGNOSTIC_ONLY_NOT_CONFIRMATION",
        "critical_caveat":"Selection uses photometric LePHARE galaxy labels, not spectroscopic/kinematic confirmation. This pair-level LOO is not the frozen full-epoch candidate estimator. No calibrated additional chromatic PSF uncertainty or revised PM significance inferred."
    }
    (out/"galaxy_only_loo_summary.json").write_text(json.dumps(result,indent=2,default=str)+"\n")
    print(json.dumps(result,indent=2,default=str))


if __name__=="__main__":
    main()
