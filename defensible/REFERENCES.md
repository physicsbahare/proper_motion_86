# References used by the defensible PM workflow

- Griggio, Nardiello & Bedin (2023), *Photometry and astrometry with JWST-II: NIRCam distortion correction*, Astronomische Nachrichten 344, e230006. DOI: 10.1002/asna.20230006; arXiv:2212.03256.
- JWST pipeline data-quality flag definitions: JWST documentation (`DO_NOT_USE`, `SATURATED`, `JUMP_DET`, `OUTLIER`).
- Photutils centroiding documentation for `centroid_2dg` and `centroid_com`.
- Astroquery MAST observation/product APIs.
- GWCS for detector-to-sky and sky-to-detector transformations using the embedded JWST CAL gWCS.

The implementation also incorporates the project's empirical validation results: exhaustive epoch search, local affine registration, Monte-Carlo centroid errors, strict DQ masking, bounded 4-pixel source association, forced-Gaussian rescue only under explicit quality cuts, and rejection of neighbor-latching PM solutions.
