"""Plot RM and polarisation maps for the Gum-Nebula SPICE-RACS DR2 subset."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.colors import LogNorm, SymLogNorm

CSV = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.rm.csv")
PSR = Path("/Users/h.qiu/cryogum/psrgum.txt")
OUT = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.rm.png")

GUM_L, GUM_B, GUM_R = 258.0, -2.0, 22.0


def load_pulsars(path: Path) -> pd.DataFrame:
    """Parse the fixed 23-token psrcat-style table; '*' -> NaN."""
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
            "name": tok[2],
            "l": f(tok[3]),
            "b": f(tok[4]),
            "ra": f(tok[5]),
            "dec": f(tok[7]),
            "dm": f(tok[11]),
            "rm": f(tok[13]),
        })
    return pd.DataFrame(rows)


df = pd.read_csv(CSV)
print(f"Loaded {len(df):,} rows")

psr = load_pulsars(PSR)
# Keep only pulsars inside the Gum cone (file is gum-targeted but be safe).
psr_sep = SkyCoord(l=psr["l"].values * u.deg, b=psr["b"].values * u.deg,
                   frame="galactic").separation(
    SkyCoord(l=GUM_L * u.deg, b=GUM_B * u.deg, frame="galactic"))
psr = psr[psr_sep.deg <= GUM_R].copy()
n_psr_rm = int(psr["rm"].notna().sum())
print(f"Pulsars in cone: {len(psr)}  (with RM: {n_psr_rm})")

# Curated RM quality cut. goodRM_flag is the SPICE-RACS-recommended flag.
good = df["goodRM_flag"] & ~df["leakage_flag"]
print(f"After goodRM & ~leakage: {good.sum():,}")
df = df[good].copy()

# Wrap longitude to [-180, 180) for symmetric plotting around l=258 (=-102).
df["l_wrap"] = ((df["l"] + 180) % 360) - 180  # not used; keep l in [0,360]

fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)

# 1) RM map ------------------------------------------------------------------
ax = axes[0, 0]
# Symmetric log-ish colour scale clipped at sensible percentile to avoid outliers.
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
ax.invert_xaxis()  # conventional: l increases to the left
# Gum cone outline
th = np.linspace(0, 2 * np.pi, 360)
ax.plot(GUM_L + GUM_R * np.cos(th) / np.cos(np.deg2rad(GUM_B)),
        GUM_B + GUM_R * np.sin(th), "k--", lw=0.8, alpha=0.6)
# Pulsars: filled stars coloured by RM where known, open stars otherwise.
with_rm = psr["rm"].notna()
ax.scatter(psr.loc[with_rm, "l"], psr.loc[with_rm, "b"],
           c=psr.loc[with_rm, "rm"], cmap="RdBu_r",
           norm=SymLogNorm(linthresh=20, vmin=-rm_clip, vmax=rm_clip, base=10),
           marker="*", s=110, edgecolor="black", linewidth=0.6, zorder=5,
           label=f"PSR w/ RM (N={int(with_rm.sum())})")
ax.scatter(psr.loc[~with_rm, "l"], psr.loc[~with_rm, "b"],
           facecolor="none", edgecolor="black",
           marker="*", s=90, linewidth=0.8, zorder=5,
           label=f"PSR no RM (N={int((~with_rm).sum())})")
ax.legend(loc="lower left", fontsize=8, framealpha=0.85)

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
ax.scatter(psr["l"], psr["b"], marker="*", s=70, facecolor="white",
           edgecolor="black", linewidth=0.6, zorder=5)

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
ax.scatter(psr["l"], psr["b"], marker="*", s=70, facecolor="white",
           edgecolor="black", linewidth=0.6, zorder=5)

# 4) RM histogram ------------------------------------------------------------
ax = axes[1, 1]
rm = df["rm"].to_numpy()
rm = rm[np.isfinite(rm)]
lo, hi = np.percentile(rm, [0.5, 99.5])
ax.hist(rm.clip(lo, hi), bins=120, color="steelblue", edgecolor="none",
        label=f"SPICE-RACS (N={len(rm):,})")
ax.axvline(0, color="k", lw=0.7)
ax.axvline(np.median(rm), color="crimson", lw=1.2,
           label=f"SPICE median = {np.median(rm):.1f} rad m$^{{-2}}$")
psr_rm = psr["rm"].dropna().to_numpy()
if psr_rm.size:
    for x in psr_rm:
        ax.axvline(np.clip(x, lo, hi), color="darkorange", lw=0.6, alpha=0.7)
    ax.axvline(np.clip(psr_rm[0], lo, hi), color="darkorange", lw=0.6,
               alpha=0.7, label=f"pulsar RM (N={psr_rm.size})")
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
