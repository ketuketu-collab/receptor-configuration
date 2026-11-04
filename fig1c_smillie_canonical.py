#!/usr/bin/env python3
"""Fig 1c dataset 2 (Smillie SCP259) — 系統別平均・z・発現マッチ帰無を1本で再生成。
2026-09-28: 旧 prep_smillie_lineage.py と監査スクリプトで系統定義が食い違っていたため統一。
確定定義（実在する23クラスタのみ）:
  B       = Follicular, Cycling B, GC
  Myeloid = Macrophages, DC1, DC2, Inflammatory Monocytes, Cycling Monocytes, CD69+/- Mast
  NK      = NKs        T = CD4*, CD8*, Tregs, ILCs, Cycling T        Plasma = Plasma
  Other   = MT-hi（ミトコンドリア高＝低品質として除外）
metadata.csv の行順は barcodes.tsv.gz と完全一致することを確認済み（assert）。"""
import gzip, os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA_ROOT
S = os.path.join(DATA_ROOT, "tlr10_uc_scrna/data/smillie_imm_export")
GENES = ["TLR10", "CD180", "LY86", "AOAH", "TLR4", "LY96", "CD14"]
LN = ["B", "Myeloid", "NK", "T", "Plasma"]
def lin(c):
    if c in ("Follicular", "Cycling B", "GC"): return "B"
    if c in ("Macrophages", "DC1", "DC2", "Inflammatory Monocytes", "Cycling Monocytes",
             "CD69+ Mast", "CD69- Mast"): return "Myeloid"
    if c == "NKs": return "NK"
    if c.startswith(("CD4", "CD8")) or c in ("Tregs", "ILCs", "Cycling T"): return "T"
    if c == "Plasma": return "Plasma"
    return "Other"
sym = pd.read_csv(f"{S}/features.tsv.gz", header=None, sep="\t")[0].astype(str).tolist()
bar = pd.read_csv(f"{S}/barcodes.tsv.gz", header=None)[0].astype(str).tolist()
meta = pd.read_csv(f"{S}/metadata.csv")
assert (meta["barcode"].astype(str).values == np.array(bar)).all(), "row order mismatch"
code = np.array([LN.index(lin(c)) if lin(c) in LN else -1 for c in meta.Cluster.astype(str)])
n = np.array([(code == i).sum() for i in range(len(LN))], dtype=float)
print("cells per lineage:", dict(zip(LN, n.astype(int))), "| Other", int((code < 0).sum()))
G = len(sym); ssum = np.zeros((G, len(LN)))
with gzip.open(f"{S}/lognorm_data.mtx.gz", "rt") as f:
    f.readline(); dims = f.readline().split()
    assert int(dims[0]) == len(bar) and int(dims[1]) == G
    for line in f:
        a = line.split(); c = int(a[0]) - 1; g = int(a[1]) - 1
        k = code[c]
        if k >= 0: ssum[g, k] += float(a[2])
M = pd.DataFrame(ssum / np.maximum(n, 1), index=sym, columns=LN)
M = M[~M.index.duplicated(keep="first")]
M.to_csv("fig1c_dataset2_lineage_means_allgenes.csv")
M.loc[[g for g in GENES if g in M.index]].to_csv("fig1c_smillie.csv")          # 7遺伝子（旧ファイルを上書き）
z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=0) + 1e-9, axis=0)
z.loc[[g for g in GENES if g in z.index]].T.to_csv("fig1c_source.csv")          # 図が読む形（lin が行）
overall = M.mean(axis=1); expressed = overall[overall > 0.005].index
rows = []
print(f"\ngenes {len(M):,} expressed {len(expressed):,}")
print(f"{'gene':7s} {'B z':>7s} {'null n':>8s} {'%ile':>8s}  verdict")
for g in ["TLR10", "CD180", "LY86", "AOAH", "TLR4", "LY96"]:
    if g not in M.index: continue
    o = overall[g]; pool = [x for x in expressed if 0.9*o <= overall[x] <= 1.1*o and x != g]
    pct = (z.loc[pool, "B"] < z.loc[g, "B"]).mean()
    v = "** specific" if pct > 0.95 else "INSIDE NULL"
    print(f"{g:7s} {z.loc[g,'B']:7.2f} {len(pool):8d} {pct:7.1%}  {v}")
    rows.append(dict(gene=g, B_z=round(z.loc[g,"B"],3), null_n=len(pool), pctile=round(pct,4), verdict=v))
pd.DataFrame(rows).to_csv("fig1c_dataset2_null.csv", index=False)
print("\nz table:\n" + z.loc[[g for g in GENES if g in z.index]].round(2).to_string())
