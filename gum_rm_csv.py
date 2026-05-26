"""Export key RM-map columns from the Gum-Nebula SPICE-RACS DR2 subset to CSV."""
from pathlib import Path

from astropy.table import Table

SRC = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.fits")
OUT = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.rm.csv")

# Minimal set for plotting an RM map + applying typical quality cuts.
COLS = [
    "cat_id",
    "ra", "dec",
    "l", "b",
    "rm", "rm_err",
    "polint", "polint_err",
    "fracpol",
    "snr_polint",
    "stokesI",
    "goodI_flag", "goodRM_flag",
    "complex_flag",
    "leakage_flag",
]

tbl = Table.read(SRC, hdu=1)
tbl[COLS].write(OUT, format="csv", overwrite=True)
print(f"Wrote {OUT} ({len(tbl):,} rows, {len(COLS)} cols)")
