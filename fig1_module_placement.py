#!/usr/bin/env python3
"""
Fig 1: human B cells swap the acyl-discrimination module.
(a) Monaco GSE107011 分取29サブセット: MD-1(LY86):MD-2(LY96) 比
(b) GSE327707 クリーン分取のベースラインCPM: B細胞はTLR4を持つがMD-2/CD14を持たない
(c) GSE214695 結腸scRNA: 3層構造（制御トリオ=B / AOAH=NK-T / TLR4-MD-2=Myeloid）
"""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import C, save, DATA_ROOT, FIGDIR

D = DATA_ROOT; OUT = FIGDIR

# ---- (a) Monaco ----
ref = pd.read_csv(f"{D}/spencer_galt/geo/GSE327417/extracted/GSM9656505_318B-2_features.tsv.gz",
                  sep="\t", header=None, names=["ens", "symbol", "type"])
m = dict(zip(ref.ens, ref.symbol))
X = pd.read_csv(f"{D}/tlr10_crosscelltype/data/GSE107011_Processed_data_TPM.txt.gz", sep="\t", index_col=0)
X.index = [str(i).split(".")[0] for i in X.index]
X["symbol"] = [m.get(i) for i in X.index]
G = ["LY86", "LY96"]
E = X[X.symbol.isin(G)].groupby("symbol").sum(numeric_only=True)
meta = pd.DataFrame({"sample": E.columns})
meta["cell"] = meta["sample"].str.split("_", n=1).str[1]
meta = meta[meta.cell != "PBMC"]
E = E[meta["sample"].values]
long = E.T.join(meta.set_index("sample")["cell"])
r = long.groupby("cell").median(numeric_only=True)
r["ratio"] = r.LY86 / r.LY96.clip(lower=1e-3)
# 両方が床の細胞種では比が無意味（NK: MD-1 1.4 / MD-2 0.1 の見かけ）。事前に floor を切る。
floor = r[["LY86", "LY96"]].max(axis=1) < 2.0   # どちらのサブユニットも2 TPM未満＝比が定義できない
print("excluded (both subunits at floor):", ", ".join(r.index[floor]))
r = r[~floor].sort_values("ratio", ascending=False)
r.to_csv(f"{OUT}/fig1a_source.csv")
print(r.round(3).to_string())

# ---- (b) GSE327707 baseline（生カウントから再計算。結果CSVには LY96 が無いため）----
GB = ["TLR4", "LY96", "CD14", "LY86", "CD180", "TLR10"]
cd = pd.read_csv(f"{D}/bcell_lipidA_whitespace/GSE327707/data/coldata.csv")
Xc = pd.read_csv(f"{D}/bcell_lipidA_whitespace/GSE327707/data/GSE327707_raw_counts.csv.gz", index_col=0)
if Xc.index.astype(str).str.startswith("ENSG").any():
    Xc.index = [m.get(str(i).split(".")[0], str(i)) for i in Xc.index]
cpm = Xc / Xc.sum() * 1e6
bl = pd.DataFrame({
    ct: cpm.loc[GB, [l for l in cd[(cd.celltype == ct) & (cd.treatment == "Baseline")].lib
                     if l in cpm.columns]].mean(axis=1)
    for ct in ["B_cell", "CD14"]})
bl.to_csv(f"{OUT}/fig1b_source.csv")
print(bl.round(2).to_string())

# ---- (c) Smillie 2019 (SCP259) lineage means ----
GC = ["TLR10", "CD180", "LY86", "AOAH", "TLR4", "LY96", "CD14"]
LN = ["B", "Myeloid", "NK", "T", "Plasma"]
lm = pd.read_csv(f"{OUT}/fig1c_smillie.csv", index_col=0).reindex(LN)[GC]
z = lm.apply(lambda c: (c - c.mean()) / (c.std(ddof=0) + 1e-9))   # 遺伝子ごとに系統間でz化
z.to_csv(f"{OUT}/fig1c_source.csv")
print(lm.round(4).to_string())
# 独立再現（GSE214695 結腸 scRNA）: 同じ順位が出るかをログに残す
pl = pd.read_csv(f"{D}/gbp2_cytosolic_arm/patient_lineage.csv")
print("\n[replication GSE214695]\n" + pl.groupby("lin")[GC].mean().round(3).to_string())

# ---- draw ----
fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.0, 1.0],
                      hspace=0.75, wspace=0.38)

ax = fig.add_subplot(gs[0, :])
sel = r.copy()
cols = []
for c in sel.index:
    if c.startswith("B_"): cols.append(C["B"])
    elif "Neutrophil" in c or "Basophil" in c: cols.append(C["gran"])
    elif "Plasmablast" in c: cols.append(C["plasma"])
    elif "monocyte" in c.lower() or c.startswith("C_mono") or "DC" in c: cols.append(C["mye"])
    else: cols.append(C["flat"])
ax.bar(range(len(sel)), sel.ratio.clip(lower=1e-3), color=cols)
ax.set_yscale("log")
ax.set_xticks(range(len(sel)))
ax.set_xticklabels(sel.index, rotation=90, fontsize=5.6)
ax.axhline(1, color="k", lw=.6, ls=":")
ax.set_ylabel("MD-1 : MD-2 (LY86/LY96, median TPM)")
ax.set_title("a  The acyl-discrimination subunit is swapped by lineage "
             "(Monaco GSE107011, sorted subsets)", loc="left")
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=c, label=l) for l, c in
                   (("B cells", C["B"]), ("myeloid/DC", C["mye"]),
                    ("granulocyte", C["gran"]), ("plasmablast", C["plasma"]),
                    ("other", C["flat"]))],
          ncol=5, loc="upper right")

ax2 = fig.add_subplot(gs[1, 0])
w = .38; xs = np.arange(len(GB))
ax2.bar(xs - w/2, bl.B_cell.values, width=w, color=C["B"], label="sorted B")
ax2.bar(xs + w/2, bl.CD14.values, width=w, color=C["mye"], label="sorted CD14 monocyte")
ax2.set_yscale("symlog", linthresh=1)
ax2.set_xticks(xs); ax2.set_xticklabels(GB, rotation=45, ha="right")
ax2.set_ylabel("baseline CPM")
ax2.set_title("b  B cells keep TLR4 but lack MD-2 and CD14\n(GSE327707, sorted, unstimulated)", loc="left")
ax2.legend(loc="upper right")

ax3 = fig.add_subplot(gs[1, 1])
im = ax3.imshow(z.values, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
ax3.set_xticks(range(len(GC))); ax3.set_xticklabels(GC, rotation=45, ha="right")
ax3.set_yticks(range(len(LN))); ax3.set_yticklabels([l.replace("_", "/") for l in LN])
ax3.set_title("c  Three tiers, three lineages\n(Smillie 2019 UC LP immune, z by gene)", loc="left")
cb = fig.colorbar(im, ax=ax3, fraction=.045, pad=.03); cb.set_label("z", fontsize=7)
save(fig, f"{OUT}/Fig1_module_placement")
