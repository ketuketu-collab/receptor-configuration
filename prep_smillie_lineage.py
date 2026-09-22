#!/usr/bin/env python3
"""Smillie 2019 (SCP259) LP免疫 54,195細胞から、配置モジュール7遺伝子の系統別平均を取り出す。
mtx を1パスで流すだけ。結果は CSV にキャッシュ（生データは read-only）。"""
import gzip, os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA_ROOT, FIGDIR

S = os.path.join(DATA_ROOT, "tlr10_uc_scrna/data/smillie_imm_export")
OUT = os.path.join(FIGDIR, "fig1c_smillie.csv")
GENES = ["TLR10", "CD180", "LY86", "AOAH", "TLR4", "LY96", "CD14"]

feat = pd.read_csv(f"{S}/features.tsv.gz", header=None, sep="\t")
sym = feat[0].astype(str).tolist()          # 1列のみ（遺伝子記号）の想定
idx = {sym[i]: i + 1 for i in range(len(sym))}   # mtx は1-origin
col2gene = {idx[g]: g for g in GENES if g in idx}
print("found:", sorted(col2gene.values()))

meta = pd.read_csv(f"{S}/metadata.csv")
def lin(c):
    if c in ("Follicular", "Cycling B", "GC"): return "B"
    if c in ("Macrophages", "DC1", "DC2", "Inflammatory Monocytes", "Cycling Monocytes"): return "Myeloid"
    if c == "NKs": return "NK"
    if c.startswith("CD4") or c.startswith("CD8") or c in ("Tregs", "ILCs", "Cycling T"): return "T"
    if c == "Plasma": return "Plasma"
    return "Other"
meta["lin"] = meta.Cluster.map(lin)
lin_arr = meta["lin"].values

acc = {g: {} for g in col2gene.values()}
with gzip.open(f"{S}/lognorm_data.mtx.gz", "rt") as fh:
    fh.readline(); fh.readline()             # header + dims
    for line in fh:
        r, c, v = line.split()
        c = int(c)
        g = col2gene.get(c)
        if g is not None:
            acc[g][int(r)] = float(v)
rows = []
for g, d in acc.items():
    vals = np.zeros(len(meta))
    if d:
        ri = np.fromiter(d.keys(), dtype=np.int64) - 1
        vv = np.fromiter(d.values(), dtype=float)
        vals[ri] = vv
    s = pd.Series(vals).groupby(lin_arr).mean()
    for k, v in s.items():
        rows.append((g, k, v))
T = pd.DataFrame(rows, columns=["gene", "lin", "mean"]).pivot(index="lin", columns="gene", values="mean")
T.to_csv(OUT)
print(T.round(4).to_string())
