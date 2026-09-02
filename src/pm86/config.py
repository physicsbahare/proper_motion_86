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

FILTER_WAVELENGTH_NM = {
    "F090W": 900,
    "F115W": 1150,
    "F150W": 1500,
    "F200W": 2000,
    "F277W": 2770,
    "F356W": 3560,
    "F410M": 4100,
    "F444W": 4440,
}

FILTER_PRIORITY = [
    "F444W",
    "F410M",
    "F356W",
    "F277W",
    "F200W",
    "F150W",
    "F115W",
    "F090W",
]
