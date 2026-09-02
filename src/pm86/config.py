from dataclasses import dataclass


@dataclass(frozen=True)
class Config:
    query_radius_arcsec: float = 3.0
    epoch_gap_days: float = 1.0
    min_pm_baseline_days: float = 30.0
    cutout_half_size: int = 200
    max_covering_exposures_per_epoch: int = 4

    target_aperture_radius_pix: float = 3.0
    target_annulus_inner_pix: float = 6.0
    target_annulus_outer_pix: float = 10.0
    target_min_snr_astrometry: float = 3.0
    target_secure_snr: float = 5.0

    control_detect_sigma: float = 4.0
    control_min_pixels: int = 6
    control_min_snr: float = 10.0
    control_target_exclusion_pix: float = 12.0
    control_edge_pix: int = 12
    min_control_matches: int = 6
    match_radius_arcsec: float = 0.35

    centroid_half_size: int = 7
    centroid_mc_draws: int = 40
    blank_apertures: int = 60

    moving_significance: float = 5.0
    zero_motion_significance: float = 3.0
    max_stationary_sigma_masyr: float = 20.0
    min_pair_direction_fraction: float = 0.75
    cross_filter_color_corr_limit: float = 0.60

    phot_hard_bad_bits: int = 1 | 2
    astrom_bad_bits: int = 1 | 2 | 4 | 16


CONFIG = Config()

# Effective/nominal wavelengths in nm. Include the NIRCam medium filters that
# commonly appear in COSMOS archival programs so they are not accidentally
# treated as wavelength=0 during pair ranking.
FILTER_WAVELENGTH_NM = {
    "F070W": 704,
    "F090W": 902,
    "F115W": 1154,
    "F140M": 1405,
    "F150W": 1501,
    "F162M": 1626,
    "F164N": 1645,
    "F182M": 1845,
    "F187N": 1874,
    "F200W": 1989,
    "F210M": 2093,
    "F212N": 2121,
    "F250M": 2503,
    "F277W": 2762,
    "F300M": 2996,
    "F322W": 3232,
    "F323N": 3237,
    "F335M": 3362,
    "F356W": 3568,
    "F360M": 3624,
    "F405N": 4052,
    "F410M": 4082,
    "F430M": 4280,
    "F444W": 4404,
    "F460M": 4624,
    "F466N": 4654,
    "F470N": 4707,
    "F480M": 4817,
}

# Deterministic tie-break order after the quantitative wavelength cost. The
# quantitative scorer, not this list, makes the main detectability decision.
FILTER_PRIORITY = [
    "F444W",
    "F430M",
    "F410M",
    "F460M",
    "F466N",
    "F470N",
    "F480M",
    "F405N",
    "F360M",
    "F356W",
    "F335M",
    "F323N",
    "F322W",
    "F300M",
    "F277W",
    "F250M",
    "F212N",
    "F210M",
    "F200W",
    "F187N",
    "F182M",
    "F164N",
    "F162M",
    "F150W",
    "F140M",
    "F115W",
    "F090W",
    "F070W",
]
