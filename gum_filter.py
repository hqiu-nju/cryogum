"""Extract SPICE-RACS DR2 sources within the Gum Nebula cone (l,b)=(258, -2), r=22 deg."""
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy import units as u

SRC = Path("/Users/h.qiu/cryogum/spice-racs.dr2.fits")
OUT = Path("/Users/h.qiu/cryogum/spice-racs.dr2.gum.fits")

GUM_CENTRE = SkyCoord(l=258 * u.deg, b=-2 * u.deg, frame="galactic")
GUM_RADIUS = 22 * u.deg

# Stream just the l, b columns first to find row indices in the cone, then
# materialise the full table for those rows only. Avoids loading 9.3M x 125 cols.
with fits.open(SRC, memmap=True) as hdul:
    data = hdul[1].data
    l = np.asarray(data["l"], dtype=np.float64)
    b = np.asarray(data["b"], dtype=np.float64)

coords = SkyCoord(l=l * u.deg, b=b * u.deg, frame="galactic")
sep = coords.separation(GUM_CENTRE)
mask = sep <= GUM_RADIUS
n_in = int(mask.sum())
print(f"Sources in Gum cone: {n_in:,} / {len(l):,} ({100 * n_in / len(l):.2f}%)")

# Now read the full table and slice. Table.read on memmap'd FITS is reasonable here.
tbl = Table.read(SRC, hdu=1)
tbl_gum = tbl[mask]
tbl_gum.meta["GUM_L"] = (258.0, "Gum cone centre Galactic longitude [deg]")
tbl_gum.meta["GUM_B"] = (-2.0, "Gum cone centre Galactic latitude [deg]")
tbl_gum.meta["GUM_R"] = (22.0, "Gum cone radius [deg]")
tbl_gum.meta["PARENT"] = (SRC.name, "Parent catalogue")

tbl_gum.write(OUT, overwrite=True)
print(f"Wrote {OUT} ({n_in:,} rows)")
