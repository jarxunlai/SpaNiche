"""Per-stage wall time of finished runs, from the files each stage writes.

Usage: python3 scripts/stage_times.py runs/A [runs/B ...]

The analysis functions save an intermediate file at the end of each stage,
so their modification times mark the stage boundaries:
  smoothing + LR matrix  -> LRintegratedmatrix.rds   (absent for celltype_and_gene)
  Moran's I              -> LRintegratedmatrix_with_moransi.rds
  NMF (init + spatial)   -> nmf_modified_fit.rds
  downstream (1b only)   -> end of the function
The first file written (spot_by_celltype_extended.csv, after < 1 s of spatial
smoothing) marks the start; start + TOTAL (from log.txt) marks the end.
"""
import os
import re
import sys


def summarise(run):
    log = open(os.path.join(run, "log.txt")).read()
    total = float(re.search(r"^TOTAL \S+: ([0-9.]+) s", log, re.M).group(1))
    start = os.path.getmtime(os.path.join(run, "spot_by_celltype_extended.csv"))
    end = start + total
    marks = [("LR matrix", "LRintegratedmatrix.rds"),
             ("Moran's I", "LRintegratedmatrix_with_moransi.rds"),
             ("NMF", "nmf_modified_fit.rds")]
    out, last = {}, start
    for name, f in marks:
        p = os.path.join(run, f)
        if os.path.exists(p):
            t = os.path.getmtime(p)
            out[name] = t - last
            last = t
    if "LR matrix" not in out:  # celltype_and_gene: smoothing + HVGs before NMF
        out = {"HVG + NMF": out.pop("NMF", 0)} | out
    out["downstream"] = max(0.0, end - last)
    return total, out, cores_in_window(run, start, end)


def cores_in_window(run, start, end):
    """Peak / mean busy cores while the analysis function itself was running
    (excludes R start-up and data loading)."""
    path = os.path.join(run, "monitor.tsv")
    rows = [l.split("\t") for l in open(path).read().splitlines()[1:]]
    t0 = os.path.getmtime(path) - float(rows[-1][0])  # monitor start
    busy = [float(c) for e, _, c in rows if start <= t0 + float(e) <= end]
    if not busy:
        return ("?", "?")
    return (f"{max(busy):.1f}", f"{sum(busy) / len(busy):.1f}")


for run in sys.argv[1:]:
    total, stages, (peak, mean) = summarise(run)
    parts = "  ".join(f"{k}={v:.0f}s" for k, v in stages.items())
    print(f"{os.path.basename(run):<22} total={total:7.1f}s  peak_cores={peak:>5} mean_cores={mean:>4}  {parts}")
