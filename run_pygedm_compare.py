import csv
import pygedm
import sys

input_file = "/Users/h.qiu/Documents/gedm_diff/gum_psrcat.csv"
output_file = "/Users/h.qiu/Documents/gedm_diff/gedm_model_comparison.csv"

rows = []
with open(input_file, "r") as f:
    reader = csv.reader(f, delimiter=";")
    headers = next(reader)  # column names
    next(reader)            # units row
    for row in reader:
        if len(row) > 2:
            rows.append(row)

results = []
skipped = []

for row in rows:
    name = row[1].strip()
    gl_str = row[9].strip()
    gb_str = row[10].strip()
    dm_str = row[19].strip()

    # Skip rows with missing coordinates or DM
    if dm_str == "*" or dm_str == "" or gl_str == "*" or gl_str == "":
        skipped.append((name, "missing DM or coordinates"))
        continue

    try:
        gl = float(gl_str)
        gb = float(gb_str)
        dm = float(dm_str)
    except ValueError:
        skipped.append((name, f"parse error: gl={gl_str}, gb={gb_str}, dm={dm_str}"))
        continue

    entry = {"NAME": name, "Gl": gl, "Gb": gb, "DM": dm}

    for method, label in [("ne2001", "NE2001"), ("ymw16", "YMW16"), ("ne2025", "NE2025")]:
        try:
            dist, tau = pygedm.dm_to_dist(gl, gb, dm, method=method)
            entry[f"dist_{label}_kpc"] = round(dist.to("kpc").value, 4)
            entry[f"tau_sc_{label}_s"]  = f"{tau.to('s').value:.6e}"
        except Exception as e:
            entry[f"dist_{label}_kpc"] = "ERROR"
            entry[f"tau_sc_{label}_s"]  = "ERROR"
            print(f"  {name} {label}: {e}", file=sys.stderr)

    # Distance discrepancies (absolute difference in kpc)
    dists = {m: entry.get(f"dist_{m}_kpc") for m in ["NE2001", "YMW16", "NE2025"]}
    if all(isinstance(v, float) for v in dists.values()):
        entry["dist_diff_NE2001_YMW16_kpc"]  = round(abs(dists["NE2001"] - dists["YMW16"]), 4)
        entry["dist_diff_NE2001_NE2025_kpc"] = round(abs(dists["NE2001"] - dists["NE2025"]), 4)
        entry["dist_diff_YMW16_NE2025_kpc"]  = round(abs(dists["YMW16"] - dists["NE2025"]), 4)
    else:
        entry["dist_diff_NE2001_YMW16_kpc"]  = "N/A"
        entry["dist_diff_NE2001_NE2025_kpc"] = "N/A"
        entry["dist_diff_YMW16_NE2025_kpc"]  = "N/A"

    # Scattering discrepancies (ratio, larger/smaller)
    def parse_tau(key):
        v = entry.get(key)
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    taus = {
        "NE2001": parse_tau("tau_sc_NE2001_s"),
        "YMW16":  parse_tau("tau_sc_YMW16_s"),
        "NE2025": parse_tau("tau_sc_NE2025_s"),
    }

    def ratio(a, b):
        if a and b and b != 0:
            r = a / b
            return round(max(r, 1/r), 4)
        return "N/A"

    entry["tau_ratio_NE2001_YMW16"]  = ratio(taus["NE2001"], taus["YMW16"])
    entry["tau_ratio_NE2001_NE2025"] = ratio(taus["NE2001"], taus["NE2025"])
    entry["tau_ratio_YMW16_NE2025"]  = ratio(taus["YMW16"],  taus["NE2025"])

    results.append(entry)

fieldnames = [
    "NAME", "Gl", "Gb", "DM",
    "dist_NE2001_kpc", "dist_YMW16_kpc", "dist_NE2025_kpc",
    "dist_diff_NE2001_YMW16_kpc", "dist_diff_NE2001_NE2025_kpc", "dist_diff_YMW16_NE2025_kpc",
    "tau_sc_NE2001_s", "tau_sc_YMW16_s", "tau_sc_NE2025_s",
    "tau_ratio_NE2001_YMW16", "tau_ratio_NE2001_NE2025", "tau_ratio_YMW16_NE2025",
]

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)

print(f"Written {len(results)} rows to {output_file}")
print(f"Skipped {len(skipped)} sources: {skipped}")
