"""RM and polarisation maps for the Gum-Nebula SPICE-RACS DR2 subset (no pulsars)."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, SymLogNorm

CSV = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.rm.csv")
OUT = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.rm.no_psr.png")

GUM_L, GUM_B, GUM_R = 258.0, -2.0, 22.0

df = pd.read_csv(CSV)
print(f"Loaded {len(df):,} rows")

good = df["goodRM_flag"] & ~df["leakage_flag"]
print(f"After goodRM & ~leakage: {good.sum():,}")
df = df[good].copy()

fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)

# 1) RM map ------------------------------------------------------------------
ax = axes[0, 0]
rm_clip = np.nanpercentile(np.abs(df["rm"]), 98)
sc = ax.scatter(
    df["l"], df["b"], c=df["rm"], s=2, alpha=0.7,
    cmap="RdBu_r",
    norm=SymLogNorm(linthresh=20, vmin=-rm_clip, vmax=rm_clip, base=10),
)
cb = fig.colorbar(sc, ax=ax, extend="both")
cb.set_label(r"RM [rad m$^{-2}$]")
ax.set_xlabel("Galactic longitude $l$ [deg]")
ax.set_ylabel("Galactic latitude $b$ [deg]")
ax.set_title(f"RM map  (N={len(df):,})")
ax.invert_xaxis()
th = np.linspace(0, 2 * np.pi, 360)
ax.plot(GUM_L + GUM_R * np.cos(th) / np.cos(np.deg2rad(GUM_B)),
        GUM_B + GUM_R * np.sin(th), "k--", lw=0.8, alpha=0.6)

# 2) Polarised intensity -----------------------------------------------------
ax = axes[0, 1]
pi = df["polint"].clip(lower=1e-6)
sc = ax.scatter(
    df["l"], df["b"], c=pi * 1e3, s=2, alpha=0.7, cmap="viridis",
    norm=LogNorm(vmin=np.nanpercentile(pi * 1e3, 2),
                 vmax=np.nanpercentile(pi * 1e3, 99)),
)
cb = fig.colorbar(sc, ax=ax, extend="both")
cb.set_label("Polarised intensity [mJy/beam]")
ax.set_xlabel("Galactic longitude $l$ [deg]")
ax.set_ylabel("Galactic latitude $b$ [deg]")
ax.set_title("Polarised intensity")
ax.invert_xaxis()

# 3) Fractional polarisation -------------------------------------------------
ax = axes[1, 0]
fp = df["fracpol"].clip(lower=1e-4, upper=1.0)
sc = ax.scatter(
    df["l"], df["b"], c=fp * 100, s=2, alpha=0.7, cmap="plasma",
    vmin=0, vmax=np.nanpercentile(fp * 100, 98),
)
cb = fig.colorbar(sc, ax=ax, extend="max")
cb.set_label("Fractional polarisation [%]")
ax.set_xlabel("Galactic longitude $l$ [deg]")
ax.set_ylabel("Galactic latitude $b$ [deg]")
ax.set_title("Fractional polarisation")
ax.invert_xaxis()

# 4) RM histogram ------------------------------------------------------------
ax = axes[1, 1]
rm = df["rm"].to_numpy()
rm = rm[np.isfinite(rm)]
lo, hi = np.percentile(rm, [0.5, 99.5])
ax.hist(rm.clip(lo, hi), bins=120, color="steelblue", edgecolor="none",
        label=f"SPICE-RACS (N={len(rm):,})")
ax.axvline(0, color="k", lw=0.7)
ax.axvline(np.median(rm), color="crimson", lw=1.2,
           label=f"median = {np.median(rm):.1f} rad m$^{{-2}}$")
ax.set_xlabel(r"RM [rad m$^{-2}$]")
ax.set_ylabel("count")
ax.set_title(f"RM distribution  (clipped to {lo:.0f}…{hi:.0f})")
ax.legend(fontsize=8)

fig.suptitle(
    f"SPICE-RACS DR2  —  Gum Nebula cone (l,b)=({GUM_L},{GUM_B}), r={GUM_R}°",
    fontsize=14,
)
fig.savefig(OUT, dpi=150)
print(f"Wrote {OUT}")
