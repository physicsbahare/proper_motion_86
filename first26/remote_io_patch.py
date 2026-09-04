"""First26-only bounded HTTP access for remote JWST CAL files.

Deep alternate-pair searches can otherwise lose an entire GitHub runner when one
MAST/fsspec range request stalls.  This patch preserves the same archive products
and FITS/gWCS path, but gives each HTTP request a finite timeout so a bad product is
recorded as an ordinary coverage-audit ERROR and the next product/pair can be tried.
"""
from __future__ import annotations

import aiohttp
from astropy.io import fits

import pm86.archive as archive

HTTP_TOTAL_TIMEOUT_S = 180
HTTP_CONNECT_TIMEOUT_S = 30
HTTP_SOCK_READ_TIMEOUT_S = 120


def _open_remote_cal_bounded(data_uri: str):
    timeout = aiohttp.ClientTimeout(
        total=HTTP_TOTAL_TIMEOUT_S,
        connect=HTTP_CONNECT_TIMEOUT_S,
        sock_connect=HTTP_CONNECT_TIMEOUT_S,
        sock_read=HTTP_SOCK_READ_TIMEOUT_S,
    )
    return fits.open(
        archive.mast_download_url(data_uri),
        mode="readonly",
        lazy_load_hdus=True,
        memmap=False,
        use_fsspec=True,
        fsspec_kwargs={
            "block_size": 4 * 1024 * 1024,
            "cache_type": "readahead",
            "client_kwargs": {"timeout": timeout},
        },
    )


def install():
    archive._open_remote_cal = _open_remote_cal_bounded
