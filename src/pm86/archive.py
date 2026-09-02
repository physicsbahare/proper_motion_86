"""MAST discovery and remote JWST CAL access.

The module deliberately discovers data from the public archive. It does not use
Betty's frozen mirror, reported astrometry, or candidate classifications.
"""

from __future__ import annotations

import copy
import re
import time
from dataclasses import dataclass
from urllib.parse import quote

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.mast import Observations
from stdatamodels import asdf_in_fits

from .config import CONFIG, FILTER_PRIORITY, FILTER_WAVELENGTH_NM


FILTER_RE = re.compile(r"F\d{3}[WMN]", re.IGNORECASE)
CAL_RE = re.compile(r"_nrc(?:a|b)(?:long|[1-4])_cal\.fits$", re.IGNORECASE)


@dataclass
class ExposureCutout:
    candidate_id: int
    ra_deg: float
    dec_deg: float
    epoch_label: str
    filter_name: str
    filename: str
    data_uri: str
    mjd: float
    detector: str
    x_full: float
    y_full: float
    x0: int
    y0: int
    sci: np.ndarray
    err: np.ndarray
    dq: np.ndarray
    gwcs: object


def normalise_filter(value: object) -> str | None:
    if value is None:
        return None
    match = FILTER_RE.search(str(value).upper())
    return match.group(0).upper() if match else None


def mast_download_url(data_uri: str) -> str:
    return (
        "https://mast.stsci.edu/api/v0.1/Download/file/?uri="
        + quote(str(data_uri), safe=":/")
    )


def _to_dataframe(table) -> pd.DataFrame:
    if table is None or len(table) == 0:
        return pd.DataFrame()
    return table.to_pandas()


def _query_region_with_retry(coord: SkyCoord, attempts: int = 5):
    """Query MAST with exponential backoff for transient connection failures."""
    last = None
    for attempt in range(attempts):
        try:
            return Observations.query_region(
                coord, radius=CONFIG.query_radius_arcsec * u.arcsec
            )
        except Exception as exc:
            last = exc
            if attempt + 1 < attempts:
                time.sleep(min(30, 2 ** attempt * 3))
    raise RuntimeError(
        f"MAST inventory query failed after {attempts} attempts: {last}"
    )


def query_candidate_inventory(candidate_id: int, ra_deg: float, dec_deg: float) -> pd.DataFrame:
    """Return all public JWST/HST imaging observations around one candidate."""
    coord = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    table = _query_region_with_retry(coord)
    df = _to_dataframe(table)

    if df.empty:
        return df

    df = df.copy()
    if "obs_collection" in df:
        df = df[df["obs_collection"].astype(str).str.upper().isin(["JWST", "HST"])]
    if "dataproduct_type" in df:
        df = df[df["dataproduct_type"].astype(str).str.lower().eq("image")]
    if "dataRights" in df:
        rights = df["dataRights"].astype(str).str.upper()
        df = df[rights.isin(["PUBLIC", "", "NAN", "NONE"])]

    df["candidate_id"] = int(candidate_id)
    df["candidate_ra_deg"] = float(ra_deg)
    df["candidate_dec_deg"] = float(dec_deg)
    df["filter_norm"] = df.get("filters", pd.Series(index=df.index, dtype=object)).map(normalise_filter)
    df["t_min_numeric"] = pd.to_numeric(df.get("t_min"), errors="coerce")
    return df.reset_index(drop=True)


def _group_epochs(obs: pd.DataFrame) -> list[dict]:
    """Group observation rows into real visits separated by > epoch_gap_days."""
    rows = obs[np.isfinite(obs["t_min_numeric"])].sort_values("t_min_numeric")
    groups: list[dict] = []
    current: list[int] = []
    previous = None

    for idx, row in rows.iterrows():
        mjd = float(row["t_min_numeric"])
        if previous is None or (mjd - previous) <= CONFIG.epoch_gap_days:
            current.append(idx)
        else:
            g = rows.loc[current]
            groups.append({"mjd": float(np.nanmedian(g["t_min_numeric"])), "rows": g.copy()})
            current = [idx]
        previous = mjd

    if current:
        g = rows.loc[current]
        groups.append({"mjd": float(np.nanmedian(g["t_min_numeric"])), "rows": g.copy()})

    return groups


def _pair_filter_cost(f1: str, f2: str) -> float:
    """Rank pair detectability for an F444W-selected candidate sample.

    A fixed preference for *any* same-filter pair can be scientifically wrong for
    very red/dropout sources: a repeated F115W epoch is methodologically neat but
    useless if the target is not detected there, while an F444W/F356W pair can
    provide real astrometry.  We therefore rank all legal pairs jointly.

    Cost is distance from F444W in wavelength for both epochs.  Same-filter pairs
    receive a modest bonus (1000 nm-equivalent), enough to prefer a nearby red
    same-filter pair but not enough for a blue nondetection pair to beat a much
    redder cross-filter alternative.  Cross-filter systematics remain explicitly
    handled downstream by the PM classifier.
    """
    w1 = FILTER_WAVELENGTH_NM.get(f1, 0)
    w2 = FILTER_WAVELENGTH_NM.get(f2, 0)
    cost = abs(w1 - 4440) + abs(w2 - 4440)
    if f1 == f2:
        cost -= 1000.0
    return float(cost)


def choose_epoch_pair(inventory: pd.DataFrame) -> dict:
    """Choose the strongest usable JWST/NIRCam epoch pair.

    All independent epoch pairs are ranked together rather than unconditionally
    preferring any same-filter pair.  This matters for the present F444W-selected
    red/dropout sample, where a blue same-filter pair can contain no measurable
    target while a red cross-filter pair does.

    Ranking:
      1. filters nearest F444W (proxy for target detectability in this sample);
      2. modest same-filter bonus to reduce chromatic/PSF systematics;
      3. longer temporal baseline as tie-breaker.

    HST rows are preserved in the inventory but are not silently mixed into the
    primary PM fit. Cross-instrument astrometry requires its own explicit test.
    """
    if inventory.empty:
        return {"status": "NO_ARCHIVE_IMAGING"}

    collection = inventory.get("obs_collection", pd.Series(index=inventory.index, dtype=object))
    jwst = inventory[collection.astype(str).str.upper().eq("JWST")].copy()
    if "instrument_name" in jwst:
        jwst = jwst[jwst["instrument_name"].astype(str).str.upper().str.contains("NIRCAM")]
    jwst = jwst[jwst["filter_norm"].notna()]

    if jwst.empty:
        return {"status": "NO_JWST_NIRCAM_IMAGING"}

    by_filter: dict[str, list[dict]] = {}
    for filt, grp in jwst.groupby("filter_norm"):
        by_filter[str(filt)] = _group_epochs(grp)

    all_epochs: list[tuple[str, dict]] = []
    for filt, epochs in by_filter.items():
        for epoch in epochs:
            all_epochs.append((filt, epoch))

    candidates = []
    for i, (f1, e1) in enumerate(all_epochs):
        for f2, e2 in all_epochs[i + 1:]:
            dt = abs(e2["mjd"] - e1["mjd"])
            if dt < CONFIG.min_pm_baseline_days:
                continue
            cost = _pair_filter_cost(f1, f2)
            # Prefer same-filter only after detectability; then prefer longer
            # baseline. FILTER_PRIORITY provides a deterministic final tie-break.
            same_penalty = 0 if f1 == f2 else 1
            p1 = FILTER_PRIORITY.index(f1) if f1 in FILTER_PRIORITY else 999
            p2 = FILTER_PRIORITY.index(f2) if f2 in FILTER_PRIORITY else 999
            candidates.append((cost, same_penalty, -dt, min(p1, p2), max(p1, p2), f1, e1, f2, e2))

    if not candidates:
        return {"status": "NO_INDEPENDENT_JWST_EPOCH_PAIR"}

    candidates.sort(key=lambda x: x[:5])
    cost, _, _, _, _, f1, e1, f2, e2 = candidates[0]

    if e1["mjd"] <= e2["mjd"]:
        early_f, early, late_f, late = f1, e1, f2, e2
    else:
        early_f, early, late_f, late = f2, e2, f1, e1

    return {
        "status": "PAIR_FOUND",
        "pair_type": "SAME_FILTER_JWST" if early_f == late_f else "CROSS_FILTER_JWST",
        "filter_early": early_f,
        "filter_late": late_f,
        "early": early,
        "late": late,
        "baseline_days_inventory": float(late["mjd"] - early["mjd"]),
        "pair_filter_cost": float(cost),
    }


def _get_products_with_retry(obs_rows: pd.DataFrame, attempts: int = 3) -> pd.DataFrame:
    last = None
    for attempt in range(attempts):
        try:
            obsids = (
                pd.to_numeric(obs_rows["obsid"], errors="coerce")
                .dropna()
                .astype("int64")
                .astype(str)
                .drop_duplicates()
                .tolist()
            )
            if not obsids:
                return pd.DataFrame()
            table = Observations.get_product_list(obsids)
            df = _to_dataframe(table)
            if df.empty:
                return df
            filename = df.get("productFilename", pd.Series(index=df.index, dtype=object)).astype(str)
            mask = filename.str.contains(CAL_RE, regex=True, na=False)
            if "productType" in df:
                mask &= df["productType"].astype(str).str.upper().eq("SCIENCE")
            df = df[mask].copy()
            if "dataURI" not in df:
                return pd.DataFrame()
            return df.drop_duplicates("dataURI").reset_index(drop=True)
        except Exception as exc:
            last = exc
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"MAST product query failed after {attempts} attempts: {last}")


def get_pair_products(pair: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    if pair.get("status") != "PAIR_FOUND":
        return pd.DataFrame(), pd.DataFrame()
    return (
        _get_products_with_retry(pair["early"]["rows"]),
        _get_products_with_retry(pair["late"]["rows"]),
    )


def _open_remote_cal(data_uri: str):
    return fits.open(
        mast_download_url(data_uri),
        mode="readonly",
        lazy_load_hdus=True,
        memmap=False,
        use_fsspec=True,
        fsspec_kwargs={"block_size": 4 * 1024 * 1024, "cache_type": "readahead"},
    )


def load_covering_cutouts(candidate_id, ra_deg, dec_deg, epoch_label, filter_name, products):
    """Read only products that truly cover the target and return local cutouts."""
    cutouts: list[ExposureCutout] = []
    audit_rows: list[dict] = []

    for _, row in products.iterrows():
        if len(cutouts) >= CONFIG.max_covering_exposures_per_epoch:
            break
        filename = str(row.get("productFilename", ""))
        data_uri = str(row.get("dataURI", ""))
        audit = {
            "candidate_id": int(candidate_id),
            "epoch_label": epoch_label,
            "requested_filter": filter_name,
            "productFilename": filename,
            "dataURI": data_uri,
            "status": "NOT_TRIED",
        }

        try:
            with _open_remote_cal(data_uri) as hdul:
                if "SCI" not in hdul or "ERR" not in hdul or "DQ" not in hdul:
                    audit["status"] = "MISSING_SCI_ERR_DQ"
                    audit_rows.append(audit)
                    continue
                ny, nx = hdul["SCI"].shape
                primary = hdul[0].header
                sci_header = hdul["SCI"].header
                with asdf_in_fits.open(hdul) as af:
                    w = copy.deepcopy(af.tree["meta"]["wcs"])
                    result = w.invert(float(ra_deg), float(dec_deg))
                    x = float(np.asarray(result[0]).squeeze())
                    y = float(np.asarray(result[1]).squeeze())

                if not (np.isfinite(x) and np.isfinite(y) and 0 <= x < nx and 0 <= y < ny):
                    audit["status"] = "OUTSIDE_DETECTOR"
                    audit_rows.append(audit)
                    continue

                half = CONFIG.cutout_half_size
                xc, yc = int(round(x)), int(round(y))
                x0, x1 = max(0, xc - half), min(nx, xc + half + 1)
                y0, y1 = max(0, yc - half), min(ny, yc + half + 1)
                sci = np.asarray(hdul["SCI"].section[y0:y1, x0:x1], dtype=np.float32)
                err = np.asarray(hdul["ERR"].section[y0:y1, x0:x1], dtype=np.float32)
                dq = np.asarray(hdul["DQ"].section[y0:y1, x0:x1], dtype=np.uint32)

                actual_filter = normalise_filter(primary.get("FILTER") or sci_header.get("FILTER")) or filter_name
                mjd = float(primary.get("EXPSTART") or sci_header.get("EXPSTART") or np.nan)
                detector = str(primary.get("DETECTOR") or sci_header.get("DETECTOR") or "UNKNOWN")

                cutouts.append(ExposureCutout(
                    candidate_id=int(candidate_id), ra_deg=float(ra_deg), dec_deg=float(dec_deg),
                    epoch_label=epoch_label, filter_name=actual_filter, filename=filename,
                    data_uri=data_uri, mjd=mjd, detector=detector, x_full=x, y_full=y,
                    x0=x0, y0=y0, sci=sci, err=err, dq=dq, gwcs=w,
                ))
                audit.update({
                    "status": "COVERED", "x_full": x, "y_full": y, "mjd": mjd,
                    "detector": detector, "actual_filter": actual_filter,
                })
                audit_rows.append(audit)
        except Exception as exc:
            audit["status"] = f"ERROR:{type(exc).__name__}:{exc}"
            audit_rows.append(audit)

    return cutouts, audit_rows
