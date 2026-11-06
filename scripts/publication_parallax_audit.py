#!/usr/bin/env python3
"""Publication audit: parallax sensitivity of two-epoch JWST motions.

IMPORTANT:
- Measures neither a new parallax nor a new independent proper motion.
- Uses an approximate EARTH barycentric observer, not JWST's exact L2
  ephemeris. Quantitative final fit requires JWST observer positions.
- Proper motion and parallax are degenerate for only two independent epochs.
- Original 2-D 'sigma' is sqrt(chi2 with 2 degrees of freedom), not a
  one-dimensional Gaussian-equivalent sigma.

Inputs:
  data/two_pm_candidates_125.csv
  data/reference_282040_centroids.csv
  optional --candidate-root path to a fresh exhaustive rerun of 1651050
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


DISTANCE_SCENARIOS_PC = {
    387132: (8.1, 19.1, 27.5),      # photometric MCMC lower/median/upper (not measured parallaxes)
    1651050: (7.2, 11.2, 17.0),     # photometric MCMC lower/median/upper (not measured parallaxes)
}
EXPECTED_PAIR = "F444W_epoch_000__F356W_epoch_000"


def tangent_basis(ra_deg: float, dec_deg: float) -> tuple[np.ndarray, np.ndarray]:
    ra, dec = np.deg2rad([ra_deg, dec_deg])
    east = np.array([-np.sin(ra), np.cos(ra), 0.0])
    north = np.array([-np.cos(ra)*np.sin(dec),
                      -np.sin(ra)*np.sin(dec), np.cos(dec)])
    return east, north


def parallax_factors_for_xyz(ra_deg, dec_deg, xyz_au):
    """Offset in tangent-plane mas for 1 mas parallax.

    Apparent offset is minus the observer's transverse barycentric vector.
    """
    e, n = tangent_basis(ra_deg, dec_deg)
    r = np.asarray(xyz_au, dtype=float)
    assert r.shape == (3,)
    return np.array([-np.dot(r, e), -np.dot(r, n)], dtype=float)


def earth_barycentric_xyz_au(mjd: float) -> np.ndarray:
    # Lazily loaded so unit tests for geometry can run independently.
    from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
    from astropy.time import Time
    import astropy.units as u
    with solar_system_ephemeris.set("builtin"):
        p = get_body_barycentric("earth", Time(mjd, format="mjd", scale="tdb"))
    return np.array([p.x.to_value(u.au), p.y.to_value(u.au),
                     p.z.to_value(u.au)])


def delta_parallax_factors(ra_deg, dec_deg, mjd_early, mjd_late):
    f1 = parallax_factors_for_xyz(
        ra_deg, dec_deg, earth_barycentric_xyz_au(mjd_early))
    f2 = parallax_factors_for_xyz(
        ra_deg, dec_deg, earth_barycentric_xyz_au(mjd_late))
    return f2-f1


def radial_significance_to_2sided_gaussian(radial_sigma: float):
    """Two-dimensional independent Gaussian null, p=exp(-r^2/2)."""
    if not np.isfinite(radial_sigma) or radial_sigma < 0:
        return np.nan, np.nan
    p = float(np.exp(-0.5*radial_sigma**2))
    return p, float(norm.isf(p/2))


def pm_after_assumed_parallax(
    mu_e, mu_n, baseline_year, delta_factors, distance_pc
):
    """Distance is assumed, not inferred. Units: mas/year, factors dimensionless."""
    pm = np.array([mu_e, mu_n], dtype=float)
    if distance_pc is None:
        return pm
    return pm - (1000.0/float(distance_pc))*np.asarray(delta_factors)/baseline_year


def epoch_times_legacy_reference(reference_csv):
    d = pd.read_csv(reference_csv)
    d["mjd"] = pd.to_numeric(d["mjd"], errors="raise")
    result = {}
    for filt in ("F444W", "F410M"):
        a = d[d["filter"].eq(filt)]["mjd"]
        if len(a) < 2:
            raise ValueError(f"Expected two or more legacy {filt} CAL exposures")
        result[filt] = float(a.mean())
    return result["F444W"], result["F410M"], len(d[d["filter"].eq("F444W")]), len(d[d["filter"].eq("F410M")])


def epoch_times_1651050(candidate_root):
    root = Path(candidate_root)
    summary = json.loads((root/"summary.json").read_text())
    if int(summary.get("candidate_id", -1)) != 1651050:
        raise ValueError("Wrong candidate root; must be candidate_1651050")
    pair_id = str(summary.get("pair_id", ""))
    if pair_id != EXPECTED_PAIR:
        raise ValueError(f"Chosen epoch pair changed to {pair_id!r}; audit before publication")
    reg_path = root/"pairs"/pair_id/"registered_target_positions.csv"
    registered = pd.read_csv(reg_path)
    # Registered table has only usable target detections.
    a = registered.loc[registered["epoch_label"].eq("early"), "mjd"]
    b = registered.loc[registered["epoch_label"].eq("late"), "mjd"]
    if len(a) < 1 or len(b) < 1:
        raise ValueError("No usable target position in both independent epochs")
    return float(a.mean()), float(b.mean()), len(a), len(b)


def build_rows(source, early_mjd, late_mjd, n_early, n_late, provenance):
    dt = (late_mjd-early_mjd)/365.25
    if abs(365.25*dt-float(source["baseline_days"])) > 2.0:
        raise ValueError(
            f"Epoch baseline disagrees with frozen result for {int(source['sample_candidate_id'])}: "
            f"{365.25*dt:.3f} vs {source['baseline_days']} days"
        )
    source_id = int(source["sample_candidate_id"])
    # Catalog coordinates are the same for both identities only after ID bookkeeping.
    coords = {
        387132: (150.1513055210011, 1.9720770310321503),
        1651050: (150.307303724545, 2.27642841317716),
    }
    ra, dec = coords[source_id]
    df = delta_parallax_factors(ra, dec, early_mjd, late_mjd)
    p_null, gaussian_equiv = radial_significance_to_2sided_gaussian(
        float(source["significance_2d"])
    )
    base = {
        "candidate_id": source_id,
        "epoch_pair": f"{source['filter_early']}/{source['filter_late']}",
        "epoch_provenance": provenance,
        "early_mjd": early_mjd,
        "late_mjd": late_mjd,
        "n_usable_early": n_early,
        "n_usable_late": n_late,
        "baseline_year": dt,
        "delta_parallax_factor_east": df[0],
        "delta_parallax_factor_north": df[1],
        "original_radial_chi2_sqrt": float(source["significance_2d"]),
        "zero_displacement_p_2d": p_null,
        "gaussian_two_sided_sigma_equivalent": gaussian_equiv,
        "sigma_mu_east_masyr": float(source["sigma_mu_alpha_cosdec_masyr"]),
        "sigma_mu_north_masyr": float(source["sigma_mu_delta_masyr"]),
        "observer": "Earth-barycentric approximation; NOT JWST L2 ephemeris",
    }
    records = []
    for kind, dist in [
        ("no_parallax", None),
        ("photometric_lower_distance", DISTANCE_SCENARIOS_PC[source_id][0]),
        ("photometric_nominal_distance", DISTANCE_SCENARIOS_PC[source_id][1]),
        ("photometric_upper_distance", DISTANCE_SCENARIOS_PC[source_id][2]),
        ("comparison_100pc", 100.0),
    ]:
        corrected = pm_after_assumed_parallax(
            float(source["mu_alpha_cosdec_masyr"]),
            float(source["mu_delta_masyr"]), dt, df, dist
        )
        records.append({
            **base, "scenario": kind,
            "assumed_distance_pc": np.nan if dist is None else float(dist),
            "assumed_parallax_mas": 0.0 if dist is None else 1000.0/dist,
            "pm_east_masyr_given_assumption": corrected[0],
            "pm_north_masyr_given_assumption": corrected[1],
            "pm_total_masyr_given_assumption": float(np.linalg.norm(corrected)),
            "parallax_displacement_east_mas": 0.0 if dist is None else 1000.0/dist*df[0],
            "parallax_displacement_north_mas": 0.0 if dist is None else 1000.0/dist*df[1],
        })
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--candidate-root", type=Path, default=None,
                    help="Fresh cached exhaustive output root/candidate_1651050")
    ap.add_argument("--output", type=Path, default=Path("publication_parallax_audit"))
    args = ap.parse_args()
    inp = args.repo_root/"data"/"two_pm_candidates_125.csv"
    sample = pd.read_csv(inp)
    if set(sample["sample_candidate_id"]) != {387132, 1651050}:
        raise ValueError("Unexpected publication candidate list")
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for _, item in sample.iterrows():
        cid = int(item["sample_candidate_id"])
        if cid == 387132:
            a, b, na, nb = epoch_times_legacy_reference(
                args.repo_root/"data"/"reference_282040_centroids.csv"
            )
            src = "independent legacy 282040 CAL exposure centroid table"
        elif args.candidate_root is not None:
            a, b, na, nb = epoch_times_1651050(args.candidate_root)
            src = "fresh exhaustive 1651050 registered target positions"
        else:
            print("SKIP candidate 1651050: fresh target-epoch MJDs required")
            continue
        records.extend(build_rows(item, a, b, na, nb, src))
    result = pd.DataFrame(records)
    if result.empty:
        raise RuntimeError("No auditable candidates")
    result.to_csv(args.output/"parallax_sensitivity.csv", index=False)
    summary = [
        "# Publication proper-motion sensitivity audit",
        "",
        "**This is not a new independent astrometric measurement.** "
        "The Earth barycentric position approximates the observer. "
        "A final joint astrometric fit must use JWST-specific ephemerides and "
        "additional independent epochs to disentangle parallax and proper motion.",
        "",
        "Two-epoch displacement is secure only relative to the adopted registration "
        "and centroiding error model. Dividing that displacement by time "
        "implicitly assumes zero differential parallax.",
        "",
    ]
    for cid, group in result.groupby("candidate_id", sort=True):
        nominal = group[group["scenario"].eq("photometric_nominal_distance")].iloc[0]
        no = group[group["scenario"].eq("no_parallax")].iloc[0]
        shift = float(np.hypot(
            nominal["parallax_displacement_east_mas"],
            nominal["parallax_displacement_north_mas"],
        ))
        summary += [
            f"## Candidate {cid}",
            f"- Epochs: MJD {no['early_mjd']:.5f} / {no['late_mjd']:.5f}, "
            f"{no['baseline_year']:.3f} years; "
            f"{int(no['n_usable_early'])}/{int(no['n_usable_late'])} usable early/late exposures.",
            f"- Published apparent annualized displacement: "
            f"({no['pm_east_masyr_given_assumption']:.2f}, "
            f"{no['pm_north_masyr_given_assumption']:.2f}) mas/yr.",
            f"- 2-D radial S/N {no['original_radial_chi2_sqrt']:.3f}; "
            f"2-D zero-displacement p = {no['zero_displacement_p_2d']:.3g}; "
            f"two-sided 1-D Gaussian equivalent = "
            f"{no['gaussian_two_sided_sigma_equivalent']:.2f} sigma "
            "(assuming independent Gaussian components, no systematics/trial correction).",
            f"- For assumed photometric distance {nominal['assumed_distance_pc']:.1f} pc, "
            f"Earth-based differential parallax shift = {shift:.1f} mas; "
            f"conditional corrected motion = "
            f"({nominal['pm_east_masyr_given_assumption']:.2f}, "
            f"{nominal['pm_north_masyr_given_assumption']:.2f}) mas/yr.",
            "- Model photometric distance is not a measured parallax. "
            "The above corrected motion is conditional on the chosen distance, "
            "not a posterior or independent verification.",
            "",
        ]
    summary += [
        "## Outstanding scientific validation",
        "- Use the JWST observer ephemeris (L2), not the Earth approximation, for final accuracy.",
        "- Obtain a third independent epoch (ideally same filter) to fit both parallax and PM.",
        "- Check cross-filter PSF/SED-related centroid bias, blends and the local registration null controls.",
        "- The F356W epoch of 1651050 has very limited usable CAL coverage; "
        "this matters independently of the formal 2-D displacement significance.",
        "- A 2-D radial statistic reported as '5 sigma' is not numerically the same "
        "as a two-sided 1-D Gaussian 5-sigma tail probability.",
        "- Do not treat other targets classified INSUFFICIENT_DATA as stationary/non-UCD.",
    ]
    (args.output/"README.md").write_text("\n".join(summary)+"\n")
    print("\n".join(summary))
    print("Wrote", args.output/"parallax_sensitivity.csv")


if __name__ == "__main__":
    main()
