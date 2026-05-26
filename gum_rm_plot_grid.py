"""Fully gridded RM / polint / fracpol maps with coloured pulsar overlays."""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
from scipy.stats import binned_statistic_2d

CSV = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.rm.csv")
PSR = Path("/Users/h.qiu/cryogum/psrgum.txt")
OUTDIR = Path("/Users/h.qiu/cryogum")

GUM_L, GUM_B, GUM_R = 258.0, -2.0, 22.0
CELL = 10 / 60    # grid cell size in degrees (10 arcmin)
MIN_PER_CELL = 2  # minimum sources per cell

PSR_SIZE = 110    # star marker size


# ---------------------------------------------------------------------------
# Helper: draw the Gum Nebula boundary circle
# ---------------------------------------------------------------------------
def gum_circle(ax):
    th = np.linspace(0, 2 * np.pi, 360)
    ax.plot(GUM_L + GUM_R * np.cos(th) / np.cos(np.deg2rad(GUM_B)),
            GUM_B + GUM_R * np.sin(th), "k--", lw=0.8, alpha=0.5)


# ---------------------------------------------------------------------------
# Helper: overlay pulsars coloured by a scalar value
# ---------------------------------------------------------------------------
def scatter_pulsars(ax, psr, col, cmap, norm, hollow_col="dimgray",
                    label_with=True, zorder=6):
    """
    Plot pulsar stars coloured by `col` (Series). NaN → hollow grey star.
    Returns the PathCollection of coloured stars (for colourbar use).
    """
    has = psr[col].notna()
    # coloured stars (with RM / flux measurement)
    sc = None
    if has.sum():
        sc = ax.scatter(
            psr.loc[has, "l"], psr.loc[has, "b"],
            c=psr.loc[has, col], cmap=cmap, norm=norm,
            marker="*", s=PSR_SIZE, edgecolor="white", linewidth=0.5,
            zorder=zorder,
            label=(f"PSR – {col} (N={has.sum()})" if label_with else None),
        )
    # hollow stars for pulsars without the measurement
    if (~has).sum():
        ax.scatter(
            psr.loc[~has, "l"], psr.loc[~has, "b"],
            facecolors="none", edgecolors=hollow_col, linewidths=0.8,
            marker="*", s=PSR_SIZE, zorder=zorder,
            label=(f"PSR – no {col} (N={(~has).sum()})" if label_with else None),
        )
    return sc


# ---------------------------------------------------------------------------
# Helper: grid a column and mask sparse cells
# ---------------------------------------------------------------------------
def make_grid(l, b, v, l_edges, b_edges, statistic="median"):
    grid, _, _, _ = binned_statistic_2d(l, b, v, statistic=statistic,
                                        bins=[l_edges, b_edges])
    cnt, _, _, _ = binned_statistic_2d(l, b, v, statistic="count",
                                       bins=[l_edges, b_edges])
    grid[cnt < MIN_PER_CELL] = np.nan
    return grid


# ---------------------------------------------------------------------------
# Load SPICE-RACS data
# ---------------------------------------------------------------------------
df = pd.read_csv(CSV)
print(f"Loaded {len(df):,} rows")
good = df["goodRM_flag"] & ~df["leakage_flag"]
df = df[good].copy()
print(f"After quality cuts: {len(df):,}")

# Grid edges (same for all panels)
l_edges = np.arange(GUM_L - GUM_R - 1, GUM_L + GUM_R + 1 + CELL, CELL)
b_edges = np.arange(GUM_B - GUM_R - 1, GUM_B + GUM_R + 1 + CELL, CELL)

grid_rm   = make_grid(df["l"], df["b"], df["rm"],              l_edges, b_edges)
grid_pi   = make_grid(df["l"], df["b"], df["polint"] * 1e3,    l_edges, b_edges)   # mJy/beam
grid_fp   = make_grid(df["l"], df["b"], df["fracpol"] * 100,   l_edges, b_edges)   # %
print(f"RM grid cells: {np.isfinite(grid_rm).sum()} / {grid_rm.size}")


# ---------------------------------------------------------------------------
# Load pulsars
# ---------------------------------------------------------------------------
def load_pulsars(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("-"):
            continue
        tok = s.split()
        if len(tok) != 23:
            continue
        def f(x):
            return np.nan if x == "*" else float(x)
        rows.append({
            "name":  tok[2],
            "l":     f(tok[3]),
            "b":     f(tok[4]),
            "rm":    f(tok[13]),   # rad m^-2
            "s1400": f(tok[19]),   # mJy  (total intensity proxy)
        })
    return pd.DataFrame(rows)


psr = load_pulsars(PSR)
psr_sep = SkyCoord(l=psr["l"].values * u.deg, b=psr["b"].values * u.deg,
                   frame="galactic").separation(
    SkyCoord(l=GUM_L * u.deg, b=GUM_B * u.deg, frame="galactic"))
psr = psr[psr_sep.deg <= GUM_R].copy().reset_index(drop=True)
print(f"Pulsars in cone: {len(psr)}  "
      f"(RM: {psr['rm'].notna().sum()}, S1400: {psr['s1400'].notna().sum()})")


# ---------------------------------------------------------------------------
# Colour norms
# ---------------------------------------------------------------------------
rm_clip = np.nanpercentile(np.abs(df["rm"]), 98)
norm_rm  = SymLogNorm(linthresh=20, vmin=-rm_clip, vmax=rm_clip, base=10)

pi_lo   = np.nanpercentile(grid_pi[np.isfinite(grid_pi)], 2)
pi_hi   = np.nanpercentile(grid_pi[np.isfinite(grid_pi)], 99)
norm_pi  = LogNorm(vmin=pi_lo, vmax=pi_hi)

fp_hi   = np.nanpercentile(grid_fp[np.isfinite(grid_fp)], 98)
norm_fp  = Normalize(vmin=0, vmax=fp_hi)

# Pulsar S1400 uses its own LogNorm (pulsar fluxes span a very different range)
s1400_vals = psr["s1400"].dropna().values
norm_s1400 = LogNorm(vmin=max(s1400_vals.min(), 0.01), vmax=s1400_vals.max()) \
    if s1400_vals.size else None


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
SUPTITLE = (f"SPICE-RACS DR2  —  Gum Nebula cone  "
            f"$(l,b)=({GUM_L},{GUM_B})$, $r={GUM_R}°$")
FIG_KW   = dict(figsize=(8, 7), constrained_layout=True)
DPI      = 150


def map_labels(ax, title):
    ax.set_xlabel("Galactic longitude $l$ [deg]")
    ax.set_ylabel("Galactic latitude $b$ [deg]")
    ax.set_title(title)
    ax.invert_xaxis()
    gum_circle(ax)


def save(fig, name):
    p = OUTDIR / name
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"Wrote {p}")


# ── 1) Gridded RM map ─────────────────────────────────────────────────────────
fig, ax = plt.subplots(**FIG_KW)
im = ax.pcolormesh(l_edges, b_edges, grid_rm.T, cmap="RdBu_r", norm=norm_rm,
                   shading="flat")
cb = fig.colorbar(im, ax=ax, extend="both")
cb.set_label(r"median RM [rad m$^{-2}$]")
map_labels(ax, f"Gridded RM  ({CELL*60:.0f}′ cells, ≥{MIN_PER_CELL}/cell)")
scatter_pulsars(ax, psr, "rm", "RdBu_r", norm_rm, label_with=True)
ax.legend(loc="lower left", fontsize=8, framealpha=0.85,
          markerscale=0.7, handlelength=1)
fig.suptitle(SUPTITLE)
save(fig, "spice-racs.dr2.gum.rm.grid.rm.png")


# ── 2) Gridded polarised intensity ────────────────────────────────────────────
fig, ax = plt.subplots(**FIG_KW)
im2 = ax.pcolormesh(l_edges, b_edges, grid_pi.T, cmap="viridis", norm=norm_pi,
                    shading="flat")
cb2 = fig.colorbar(im2, ax=ax, extend="both")
cb2.set_label("median pol. intensity [mJy/beam]")
map_labels(ax, "Gridded polarised intensity")
scatter_pulsars(ax, psr, "s1400", "viridis", norm_pi, label_with=True)
ax.legend(loc="lower left", fontsize=8, framealpha=0.85,
          markerscale=0.7, handlelength=1)
fig.suptitle(SUPTITLE)
save(fig, "spice-racs.dr2.gum.rm.grid.polint.png")


# ── 3) Gridded fractional polarisation ────────────────────────────────────────
fig, ax = plt.subplots(**FIG_KW)
im3 = ax.pcolormesh(l_edges, b_edges, grid_fp.T, cmap="plasma", norm=norm_fp,
                    shading="flat")
cb3 = fig.colorbar(im3, ax=ax, extend="max")
cb3.set_label("median fractional pol. [%]")
map_labels(ax, "Gridded fractional polarisation")
scatter_pulsars(ax, psr, "rm", "RdBu_r", norm_rm, label_with=True)
ax.legend(loc="lower left", fontsize=8, framealpha=0.85,
          markerscale=0.7, handlelength=1)
fig.suptitle(SUPTITLE)
save(fig, "spice-racs.dr2.gum.rm.grid.fracpol.png")


# ── 4) RM histogram ───────────────────────────────────────────────────────────
from matplotlib.lines import Line2D
fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
rm_arr = df["rm"].to_numpy()
rm_arr = rm_arr[np.isfinite(rm_arr)]
lo, hi = np.percentile(rm_arr, [0.5, 99.5])
ax.hist(rm_arr.clip(lo, hi), bins=120, color="steelblue", edgecolor="none")
ax.axvline(0, color="k", lw=0.7)
ax.axvline(np.median(rm_arr), color="crimson", lw=1.2)
psr_rm_vals = psr["rm"].dropna().to_numpy()
if psr_rm_vals.size:
    cmap_psr = plt.get_cmap("RdBu_r")
    for x in psr_rm_vals:
        ax.axvline(np.clip(x, lo, hi), color=cmap_psr(norm_rm(x)),
                   lw=0.9, alpha=0.75, zorder=3)
ax.legend(handles=[
    Line2D([0], [0], color="steelblue", lw=6),
    Line2D([0], [0], color="crimson", lw=1.5),
    Line2D([0], [0], color="gray", lw=1.2, alpha=0.75),
], labels=[
    f"SPICE-RACS (N={len(rm_arr):,})",
    f"SPICE median = {np.median(rm_arr):.1f} rad m$^{{-2}}$",
    f"pulsar RM (N={psr_rm_vals.size}, coloured by RM)",
], fontsize=9, loc="upper left")
ax.set_xlabel(r"RM [rad m$^{-2}$]")
ax.set_ylabel("count")
ax.set_title(f"RM distribution  (clipped to {lo:.0f} … {hi:.0f})")
fig.suptitle(SUPTITLE)
save(fig, "spice-racs.dr2.gum.rm.grid.hist.png")
