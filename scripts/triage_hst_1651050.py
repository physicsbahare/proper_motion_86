#!/usr/bin/env python3
"""Conservative HST 1651050 archive-duplicate and dither repeatability QA.

Reads the existing original-pixel screening CSV. Rootnames are *exposures*,
not independent astrometric epochs. No astrometric centroid is validated here.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd


def inferred_filter(row):
    # WFPC2 primary FITS FILTER can be an integer code, e.g. "10".
    f = str(row.get("filter", "")).upper()
    if re.fullmatch(r"F\d{3,4}[WMN]", f):
        return f
    fn = str(row.get("filename", "")).upper()
    m = re.search(r"F\d{3,4}[WMN]", fn)
    if m:
        return m.group(0)
    if f.startswith("G"):
        return f + " (GRISM)"
    return "UNKNOWN"  # do not invent the physical filter from a numeric code


def triage(df):
    covered = df.loc[df["coverage_status"].eq("COVERS_WITH_STAMP")].copy()
    if covered.empty:
        return pd.DataFrame(), pd.DataFrame()
    covered["snr"] = pd.to_numeric(covered.get("conservative_snr"), errors="coerce")
    covered["bad"] = pd.to_numeric(covered.get("bad_aperture_pixels"), errors="coerce")
    covered["chi2"] = pd.to_numeric(covered.get("gauss_chi2_reduced"), errors="coerce")
    covered["sigma"] = pd.to_numeric(covered.get("gauss_sigma_pix"), errors="coerce")
    covered["inferred_filter"] = covered.apply(inferred_filter, axis=1)
    rows=[]
    for root, g in covered.groupby("rootname", dropna=False):
        valid=g["snr"].dropna()
        # Prefer explicit filter names in an alias over an uninformative primary code.
        names=[f for f in g.inferred_filter if f!="UNKNOWN"]
        filter_name=names[0] if names else "UNKNOWN"
        high=bool((valid>=5).any())
        dq_bad=bool((g.bad.fillna(0)>0).any())
        any_failed=bool(g.snr.isna().any())
        sigma_floor=bool((g.sigma.dropna()<=0.700001).any())
        # This is deliberately a *screening* classification, never astrometric approval.
        if high and (dq_bad or any_failed or sigma_floor):
            status="UNVALIDATED_HIGH_SNR_DQ_OR_FIT_ANOMALY"
        elif high:
            status="HIGH_SNR_REQUIRES_PSF_AND_REGISTRATION"
        else:
            status="NO_RELIABLE_DETECTION_IN_THIS_ROOT"
        rows.append(dict(rootname=root, filter_name=filter_name,
            expstart_mjd=float(pd.to_numeric(g.expstart_mjd,errors="coerce").median()),
            archive_alias_rows=len(g), valid_forced_measurements=len(valid),
            forced_snr_min=float(valid.min()) if len(valid) else np.nan,
            forced_snr_max=float(valid.max()) if len(valid) else np.nan,
            max_bad_aperture_pixels=float(g.bad.max()) if g.bad.notna().any() else np.nan,
            max_gaussian_reduced_chi2=float(g.chi2.max()) if g.chi2.notna().any() else np.nan,
            gaussian_sigma_at_lower_bound=sigma_floor,
            inconsistent_alias_measurement=any_failed or (len(valid)>1 and (valid.max()-valid.min()>0.5)),
            screening_status=status))
    roots=pd.DataFrame(rows).sort_values(["expstart_mjd","rootname"])
    bands=[]
    for filt,g in roots.groupby("filter_name"):
        clean=g[(g.forced_snr_max>=5)&(g.max_bad_aperture_pixels.fillna(1)==0)
            &~g.gaussian_sigma_at_lower_bound&~g.inconsistent_alias_measurement]
        bands.append(dict(filter_name=filt, unique_exposure_roots=len(g),
            roots_with_any_forced_snr_ge_5=int((g.forced_snr_max>=5).sum()),
            roots_with_clean_repeatable_snr_ge_5=len(clean),
            repeatability_status="NO_VALIDATED_ASTROMETRIC_DETECTION"))
    return roots, pd.DataFrame(bands).sort_values("filter_name")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.input)
    roots,bands=triage(df)
    roots.to_csv(out/"HST_UNIQUE_EXPOSURE_TRIAGE.csv",index=False)
    bands.to_csv(out/"HST_FILTER_REPEATABILITY.csv",index=False)
    lines=["# HST archive aliases and independent-root screening for 1651050","",
        f"WCS-covering product rows: {int(df.coverage_status.eq('COVERS_WITH_STAMP').sum())}.",
        f"Unique HST exposure rootnames: {len(roots)}.",
        "A repeated rootname is a duplicate archive/calibration product, not an independent exposure.",
        "Even independent exposures within one visit do not establish independent parallax phases.","",
        "## Per-filter repeatability (NOT astrometry)","",
        bands.to_markdown(index=False) if len(bands) else "No covering exposures.","",
        "## Per-root diagnostics", "",
        roots.to_markdown(index=False) if len(roots) else "None.","",
        "Any high forced-aperture S/N with DQ contamination, a Gaussian width at its",
        "hard lower bound, or non-reproduction across dithers remains an artifact suspect.",
        "No third independent registered centroid, intrinsic PM or joint parallax is established."]
    (out/"HST_ANOMALY_TRIAGE.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines),flush=True)


if __name__=="__main__":main()
