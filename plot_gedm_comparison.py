import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.skyview import SkyView
from reproject import reproject_interp

input_file  = "/Users/h.qiu/Documents/gedm_diff/gedm_model_comparison.csv"
output_file = "/Users/h.qiu/Documents/gedm_diff/gedm_dist_diff_YMW16_NE2025.png"

df = pd.read_csv(input_file)
df = df[pd.to_numeric(df["dist_diff_YMW16_NE2025_kpc"], errors="coerce").notna()].copy()
df["dist_diff_YMW16_NE2025_kpc"] = df["dist_diff_YMW16_NE2025_kpc"].astype(float)

# coordinate extent with a small margin
gl_min = df["Gl"].min() - 3
gl_max = df["Gl"].max() + 3
gb_min = df["Gb"].min() - 3
gb_max = df["Gb"].max() + 3
gl_cen = (gl_min + gl_max) / 2
gb_cen = (gb_min + gb_max) / 2

# --- fetch H-alpha image in native RA/Dec WCS ---
print("Fetching SHASSA H-alpha image from SkyView...")
centre_eq = SkyCoord(l=gl_cen * u.deg, b=gb_cen * u.deg, frame="galactic").icrs
fits_list = SkyView.get_images(
    position=centre_eq,
    survey=["SHASSA H"],
    width=(gl_max - gl_min) * u.deg,
    height=(gb_max - gb_min) * u.deg,
    pixels="900,600",
)
hdu = fits_list[0][0]
print(f"Fetched image shape: {hdu.data.shape}")

# mask the integer BLANK sentinel (-2^31)
hdu_data = hdu.data.astype(np.float32)
hdu_data[hdu_data < -1e9] = np.nan
hdu.data = hdu_data

# --- target galactic WCS (CAR projection) ---
nx, ny = 900, 600
target_wcs = WCS(naxis=2)
target_wcs.wcs.crpix = [nx / 2 + 0.5, ny / 2 + 0.5]
target_wcs.wcs.cdelt = [(gl_max - gl_min) / nx, (gb_max - gb_min) / ny]
target_wcs.wcs.crval = [gl_cen, gb_cen]
target_wcs.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]

print("Reprojecting to galactic coordinates...")
halpha_reproj, _ = reproject_interp(hdu, target_wcs, shape_out=(ny, nx))

# asinh stretch: compresses the dynamic range while preserving faint structure
a = np.nanpercentile(halpha_reproj, 1)
b = np.nanpercentile(halpha_reproj, 99.5)
halpha_display = np.arcsinh((halpha_reproj - a) / (b - a) * 10) / np.arcsinh(10)
halpha_display = np.clip(halpha_display, 0, 1)

# --- plot ---
fig, ax = plt.subplots(figsize=(13, 5.5))

ax.imshow(
    halpha_display,
    origin="lower",
    extent=[gl_max, gl_min, gb_min, gb_max],
    cmap="Blues_r",
    vmin=0, vmax=1,
    aspect="auto",
    interpolation="bilinear",
    zorder=1,
)

norm = mcolors.LogNorm(
    vmin=df["dist_diff_YMW16_NE2025_kpc"].replace(0, np.nan).min(),
    vmax=df["dist_diff_YMW16_NE2025_kpc"].max(),
)

sc = ax.scatter(
    df["Gl"], df["Gb"],
    c=df["dist_diff_YMW16_NE2025_kpc"],
    cmap="cividis",   # perceptually uniform, safe for all colour-vision deficiencies
    norm=norm,
    s=45,
    edgecolors="white",
    linewidths=0.4,
    zorder=3,
)

cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label(r"$|\Delta d|$ YMW16 $-$ NE2025 (kpc)", fontsize=11)

ax.set_xlim(gl_max, gl_min)
ax.set_ylim(gb_min, gb_max)
ax.set_xlabel("Galactic Longitude $l$ (deg)", fontsize=12)
ax.set_ylabel("Galactic Latitude $b$ (deg)", fontsize=12)
ax.set_title(
    r"Distance discrepancy: YMW16 vs NE2025  |  H$\alpha$ background (SHASSA, SkyView)",
    fontsize=12,
)
ax.grid(True, linestyle="--", alpha=0.3, color="white")

fig.tight_layout()
fig.savefig(output_file, dpi=150)
print(f"Saved to {output_file}")
