#!/usr/bin/env python3
"""JWST-L2 (JPL Horizons -170) conditional differential-parallax audit.

This is a **sensitivity test**: two astrometric epochs do NOT independently
solve for parallax and two-component proper motion simultaneously.
The nominal, lower, and upper photometric distances are assumed, not fitted.
Horizons vectors must be JWST SSB ICRF/Earth-equatorial vectors; failures must
not silently fall back to Earth's barycentric position.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from astropy.time import Time
from astroquery.jplhorizons import Horizons

RA=150.307303724545
DEC=2.27642841317716
DISTANCES_PC=(7.2,11.2,17.0)
JWST_ID="-170"


def basis_icrs(ra_deg=RA,dec_deg=DEC):
    ra,dec=np.deg2rad([ra_deg,dec_deg])
    east=np.array([-np.sin(ra),np.cos(ra),0.])
    north=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
    return east,north


def vectors_to_parallax_factors(xyz_au,ra_deg=RA,dec_deg=DEC):
    """A 1-mas parallax gives offset -R_perp mas where R in AU."""
    xyz=np.asarray(xyz_au,dtype=float)
    east,north=basis_icrs(ra_deg,dec_deg)
    return -np.column_stack([xyz@east,xyz@north])


def get_jwst_vectors(mjd):
    """Solar-system barycentric JWST observer coordinates (AU), ICRF axes.

    astroquery uses refplane='earth' for equatorial reference, per Horizons
    documentation. Explicit '-170' selects JWST, *not* Earth's center.
    """
    dates=Time(np.asarray(mjd,dtype=float),format="mjd",scale="utc").tdb.jd
    q=Horizons(id=JWST_ID,location="@0",epochs=list(dates))
    tab=q.vectors(refplane="earth",aberrations="geometric")
    xyz=np.array([[float(tab["x"][i]),float(tab["y"][i]),float(tab["z"][i])]
                  for i in range(len(tab))])
    if xyz.shape!=(len(mjd),3) or not np.isfinite(xyz).all():
        raise ValueError(f"Wrong Horizons vector shape: {xyz.shape}")
    if not np.all((np.linalg.norm(xyz,axis=1)>.7)&(np.linalg.norm(xyz,axis=1)<1.3)):
        raise ValueError("JWST barycentric norm unexpected. Check Horizons coordinate origin.")
    return xyz,{"target":JWST_ID,"observer":"@0 Solar-System barycenter",
                "refplane":"earth = equatorial ICRF",
                "aberration":"geometric","epochs_TDB_JD":list(map(float,dates))}


def conditional_pm(original_pm,err_pm,delta_f,baseline_year,distance_pc):
    parallax=1000./float(distance_pc)
    correction=parallax*np.asarray(delta_f)/baseline_year
    residual=np.asarray(original_pm,dtype=float)-correction
    err=np.asarray(err_pm,dtype=float)
    radial=np.sqrt(np.sum((residual/err)**2))
    p=float(chi2.sf(radial**2,df=2))
    z=float(norm.isf(p/2)) if p>0 else float("inf")
    return dict(assumed_distance_pc=distance_pc,assumed_parallax_mas=parallax,
                parallax_displacement_east_mas=parallax*delta_f[0],
                parallax_displacement_north_mas=parallax*delta_f[1],
                conditional_pm_east_masyr=residual[0],
                conditional_pm_north_masyr=residual[1],
                residual_2d_sqrt_chi2=radial,
                residual_p_2d_zero_motion=p,
                residual_gaussian_1d_two_sided_equivalent=z)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--epochs",default="data/1651050_astrometric_epoch_manifest.csv")
    ap.add_argument("--summary",default="data/two_pm_candidates_125.csv")
    ap.add_argument("--output",default="publication_check/jwst_parallax")
    args=ap.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    epochs=pd.read_csv(args.epochs)
    filt=epochs.loc[epochs["astrometrically_usable"].astype(bool)].sort_values("mjd")
    if len(filt)!=3 or filt["filter"].tolist()!=["F444W","F444W","F356W"]:
        raise ValueError("Unexpected 1651050 accepted exposure manifest")
    try: xyz,provenance=get_jwst_vectors(filt["mjd"].to_numpy())
    except Exception as exc:
        (out/"HORIZONS_ERROR.txt").write_text(str(exc)+"\nNo Earth fallback was applied.\n")
        raise
    factor=vectors_to_parallax_factors(xyz)
    f1=factor[:2].mean(axis=0);f2=factor[2];delta=f2-f1
    mjd1=float(filt["mjd"][:2].mean());mjd2=float(filt["mjd"].iloc[2])
    dt=(mjd2-mjd1)/365.25
    src=pd.read_csv(args.summary)
    row=src.loc[src.sample_candidate_id.eq(1651050)].iloc[0]
    if abs((mjd2-mjd1)-float(row.baseline_days))>2:
        raise RuntimeError("Rerun and archived astrometric baselines disagree")
    mu=np.array([row.mu_alpha_cosdec_masyr,row.mu_delta_masyr],dtype=float)
    err=np.array([row.sigma_mu_alpha_cosdec_masyr,row.sigma_mu_delta_masyr],dtype=float)
    records=[]
    for distance in DISTANCES_PC:
        records.append(conditional_pm(mu,err,delta,dt,distance))
    if not np.isfinite(delta).all():
        raise ValueError("Non-finite JWST differential parallax factor")
    df=pd.DataFrame(records)
    df.to_csv(out/"jwst_l2_parallax_conditional_pm.csv",index=False)
    rec=filt[["filename","filter","mjd"]].copy()
    rec[["observer_x_au","observer_y_au","observer_z_au"]]=xyz
    rec[["factor_east","factor_north"]]=factor
    rec.to_csv(out/"jwst_horizons_observer_vectors.csv",index=False)
    receipt=dict(provenance,ra_deg=RA,dec_deg=DEC,mjd_early=mjd1,mjd_late=mjd2,
        baseline_days=mjd2-mjd1,delta_parallax_factor_east=float(delta[0]),
        delta_parallax_factor_north=float(delta[1]),
        apparent_pm_east_masyr=float(mu[0]),apparent_pm_north_masyr=float(mu[1]),
        formal_sigma_pm_east_masyr=float(err[0]),formal_sigma_pm_north_masyr=float(err[1]),
        warning="Conditional 2-epoch solution; distance fixed, errors assume independent axes; no calibration systematics.")
    (out/"PROVENANCE.json").write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))
    print(df.to_string(index=False))
    summary=[
        "# Candidate 1651050: JWST-L2 differential parallax sensitivity",
        "",
        "JPL Horizons target **-170** was queried relative to solar-system barycenter "
        "with equatorial ICRF axes. No approximation using Earth coordinates is used.",
        "",
        f"Early F444W MJD {mjd1:.6f}, late F356W MJD {mjd2:.6f}; "
        f"baseline {(mjd2-mjd1):.4f} days.",
        "",
        f"Differential parallax factors east={delta[0]:.6f}, north={delta[1]:.6f}.",
        "",
        "The table contains *conditional* proper motion after subtracting the "
        "parallax predicted by an assumed photometric distance (7.2, 11.2, 17 pc). "
        "It is not a measured trigonometric parallax.",
        "",
        "Two epochs are not sufficient to independently determine parallax and "
        "proper motion; correlated astrometric systematics are also not in the formal "
        "2-D Gaussian p-values. A 2-D radial statistic must not be called a "
        "one-dimensional Gaussian n-sigma value.",
        "",
        "Cross-filter PSF and blending must be assessed independently before "
        "assigning a secure-motion classification.",
    ]
    (out/"RESULTS_README.md").write_text("\n".join(summary)+"\n")


if __name__=="__main__":
    main()
