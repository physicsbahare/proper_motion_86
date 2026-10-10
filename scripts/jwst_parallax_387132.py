#!/usr/bin/env python3
"""JWST L2 conditional parallax check for validated candidate 387132/legacy 282040.

Source: archived six original CAL centroids from data/reference_282040_centroids.csv.
No new positions, proper motion, or independent astrometric parallax are measured.
"""
import json
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scripts.jwst_parallax_1651050 import (
    get_jwst_vectors, vectors_to_parallax_factors,
    conditional_pm, parallax_only_fit_2d,
)

RA=150.1513055210011
DEC=1.9720770310321503


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",default="publication_check/jwst_parallax_387132")
    args=ap.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    obs=pd.read_csv("data/reference_282040_centroids.csv")
    if len(obs)!=6 or obs["filter"].value_counts().to_dict()!={"F444W":4,"F410M":2}:
        raise ValueError("Unexpected validated 6-CAL exposure reference table")
    xyz,provenance=get_jwst_vectors(obs["mjd"].to_numpy(float))
    factors=vectors_to_parallax_factors(xyz,RA,DEC)
    obs[["jwst_x_au","jwst_y_au","jwst_z_au"]]=xyz
    obs[["parallax_factor_east","parallax_factor_north"]]=factors
    obs.to_csv(out/"jwst_horizons_vectors_six_exposures.csv",index=False)
    e=obs["filter"].eq("F444W")
    l=obs["filter"].eq("F410M")
    dfactors=factors[l].mean(axis=0)-factors[e].mean(axis=0)
    dt=(obs.loc[l,"mjd"].mean()-obs.loc[e,"mjd"].mean())/365.25
    frozen=pd.read_csv("data/two_pm_candidates_125.csv")
    row=frozen.loc[frozen.sample_candidate_id.eq(387132)].iloc[0]
    if abs(float(row.baseline_days)-dt*365.25)>1:
        raise RuntimeError("Reference baseline mismatch")
    mu=np.array([row.mu_alpha_cosdec_masyr,row.mu_delta_masyr],float)
    sd=np.array([row.sigma_mu_alpha_cosdec_masyr,row.sigma_mu_delta_masyr],float)
    rows=[conditional_pm(mu,sd,dfactors,dt,dist) for dist in (8.1,19.1,27.5)]
    pd.DataFrame(rows).to_csv(out/"387132_jwst_l2_conditional_motion.csv",index=False)
    zero=parallax_only_fit_2d(mu,sd,dfactors,dt)
    receipt={"source":"6 archived CAL centroids, independent already-validated 282040 PM solution",
             "candidate_id":387132,"legacy_id":282040,"ra_deg":RA,"dec_deg":DEC,
             "differential_parallax_factors":list(map(float,dfactors)),
             "baseline_years":dt,"horizons":provenance,
             "zero_pm_diagnostic":zero,
             "warning":"Photometric distances are assumed. Two epochs do not solve parallax and PM independently."}
    (out/"PROVENANCE_AND_QUALITY.json").write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__=="__main__":
    main()
