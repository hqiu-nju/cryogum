# cryogum

Preparatory analysis for CryoPAF proposals, focusing on the polarised radio sky in the direction of the **Gum Nebula** using SPICE-RACS DR2 data.

## Overview

The Gum Nebula (centred at *(l, b)* = (258°, −2°), radius ≈ 22°) is a large, nearby H II region whose ionised gas is expected to imprint a distinctive Faraday rotation signature on background polarised sources. This repository contains scripts to:

- Extract the SPICE-RACS DR2 source catalogue within the Gum Nebula cone
- Compute and grid rotation measure (RM), polarised intensity, and fractional polarisation maps at 10′ resolution
- Overlay ATNF pulsar catalogue entries with colour-coded RM and flux measurements
- Compare DM-derived distances from the YMW16 and NE2001 electron density models

## Repository structure

```
cryogum/
├── gum_rm_csv.py            # Extract Gum Nebula cone subset from SPICE-RACS DR2 FITS
├── gum_filter.py            # Quality-cut helpers (goodRM flag, leakage flag)
├── gum_rm_plot.py           # 4-panel scatter plot: RM, polint, fracpol, histogram
├── gum_rm_plot_no_psr.py    # Same scatter plot without pulsar overlays
├── gum_rm_plot_grid.py      # Gridded maps (10′ cells) — 4 separate output figures
├── run_pygedm_compare.py    # Run pygedm for YMW16 and NE2001 on cone pulsars
├── plot_gedm_comparison.py  # Plot YMW16 vs NE2001 DM distance comparison
├── psrgum.txt               # ATNF pulsar catalogue extract for the Gum Nebula cone
└── *.png                    # Output plots (see below)
```

## Output plots

| File | Description |
|---|---|
| `spice-racs.dr2.gum.rm.png` | 4-panel scatter map (RM, polint, fracpol, histogram) |
| `spice-racs.dr2.gum.rm.no_psr.png` | Same, without pulsars |
| `spice-racs.dr2.gum.rm.grid.rm.png` | Gridded RM map with pulsar overlays |
| `spice-racs.dr2.gum.rm.grid.polint.png` | Gridded polarised intensity map |
| `spice-racs.dr2.gum.rm.grid.fracpol.png` | Gridded fractional polarisation map |
| `spice-racs.dr2.gum.rm.grid.hist.png` | RM distribution histogram |
| `gedm_dist_diff_YMW16_NE2025.png` | YMW16 vs NE2001 DM distance comparison |

## Data requirements

The following large data files are required locally but are excluded from the repository (>100 MB):

| File | Description |
|---|---|
| `spice-racs.dr2.fits` | Full SPICE-RACS DR2 source catalogue |
| `spice-racs.dr2.gum.fits` | Gum Nebula cone subset (FITS) |
| `spice-racs.dr2.gum.rm.csv` | Gum Nebula cone subset with RM columns (CSV) |
| `RACS-low3_INITIAL_*.fits` | RACS-low3 component/source catalogues |

## Environment

Scripts require the `astro_fit` conda environment:

```bash
conda activate astro_fit
python gum_rm_plot_grid.py
```

Key dependencies: `numpy`, `pandas`, `matplotlib`, `astropy`, `scipy`, `pygedm`

## Context

This work supports preparatory science for the **CryoPAF** (Cryogenic Phased Array Feed) instrument proposals, investigating the Faraday rotation structure of the Gum Nebula as a foreground screen.
