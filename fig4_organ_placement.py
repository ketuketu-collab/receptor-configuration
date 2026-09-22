#!/usr/bin/env python3
"""
Fig 4: configuration varies by organ and by cell type.
(a) TICA 54,934 B細胞: 臓器別 TLR10（全身性 > 腸管）。サブタイプ固定でも同方向。
(b) HPA 154細胞種: 上皮にはモジュールが無い（CASP4だけある）。
"""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import C, save, DATA_ROOT, FIGDIR
from scipy.stats import mannwhitneyu
import scanpy as sc

D = DATA_ROOT; OUT = FIGDIR
# 事前定義（記憶の分類に合わせる）: 腸管系8臓器（腸間膜リンパ節を含む） vs 全身性リンパ組織5臓器
GUT = {"DUO", "JEJEPI", "JEJLP", "ILE", "CAE", "TCL", "SCL", "MLN"}
SYS = {"SPL", "LLN", "THY", "BMA", "BLD"}

A = sc.read_h5ad(f"{D}/tica_pan_immune_bcell/data/b-cells.h5ad", backed="r")
obs = A.obs.copy()
tis_col = [c for c in ("Tissue", "tissue", "Organ", "organ") if c in obs.columns][0]
ct_col = "Manually_curated_celltype"
gi = list(A.var_names).index("TLR10")
x = A[:, gi].to_memory().X
x = np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
obs["TLR10"] = x
obs["grp"] = np.where(obs[tis_col].astype(str).isin(GUT), "gut",
              np.where(obs[tis_col].astype(str).isin(SYS), "systemic", "other"))
obs.to_csv(f"{OUT}/fig4a_cells.csv.gz", columns=[tis_col, ct_col, "TLR10", "grp"], index=False)

by_t = obs.groupby(tis_col).TLR10.agg(["mean", lambda v: (v > 0).mean() * 100, "size"])
by_t.columns = ["mean", "pct_pos", "n"]
by_t = by_t[by_t.n >= 100].sort_values("mean", ascending=False)
print(by_t.round(3).to_string())

res = []
for sub in ["Memory B cells", "Naive B cells"]:
    o = obs[obs[ct_col].astype(str) == sub]
    if len(o) == 0: continue
    g = o[o.grp == "gut"].TLR10; s = o[o.grp == "systemic"].TLR10
    res.append((sub, g.mean(), (g > 0).mean() * 100, len(g),
                s.mean(), (s > 0).mean() * 100, len(s), mannwhitneyu(g, s).pvalue))
R = pd.DataFrame(res, columns=["subset", "gut_mean", "gut_pct", "gut_n",
                               "sys_mean", "sys_pct", "sys_n", "p"])
R.to_csv(f"{OUT}/fig4a_subset.csv", index=False)
print(R.to_string())
print("subset labels:", obs[ct_col].value_counts().head(12).to_string())

# ---- (b) HPA ----
H = pd.read_csv(f"{D}/md1_ladder/rna_single_cell_type.tsv", sep="\t")
GB = ["TLR4", "LY96", "LY86", "CD180", "CASP4"]
h = H[H["Gene name"].isin(GB)].pivot_table(index="Cell type", columns="Gene name", values="nCPM")
CT = ["paneth cells", "tuft cells", "enterocytes", "mucous neck cells",
      "b-cells", "plasma cells", "macrophages", "monocytes", "neutrophils"]
ct_map = {c.lower(): c for c in h.index}
sel = [ct_map[c] for c in CT if c in ct_map]
print("missing cell types:", [c for c in CT if c not in ct_map])
hb = h.loc[sel, GB]
hb.to_csv(f"{OUT}/fig4b_source.csv")
print(hb.round(1).to_string())

# ---- draw ----
fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.45, 1.0],
                      hspace=0.85, wspace=0.35)

ax = fig.add_subplot(gs[0, :])
cols = [C["B"] if t in SYS else (C["epi"] if t in GUT else C["flat"]) for t in by_t.index]
ax.bar(range(len(by_t)), by_t["pct_pos"], color=cols)
ax.set_xticks(range(len(by_t))); ax.set_xticklabels(by_t.index, rotation=90, fontsize=6)
ax.set_ylabel("TLR10$^+$ B cells (%)")
ax.set_title("a  TLR10 is highest in systemic lymphoid tissue, lowest in gut "
             "(TICA, 54,934 B cells, 17 tissues)", loc="left")
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=c, label=l) for l, c in
                   (("systemic lymphoid", C["B"]), ("gut-associated (incl. MLN)", C["epi"]),
                    ("other", C["flat"]))], ncol=3, loc="upper right")

ax2 = fig.add_subplot(gs[1, 0])
w = .35; xs = np.arange(len(R))
ax2.bar(xs - w/2, R.gut_pct, width=w, color=C["epi"], label="gut-associated")
ax2.bar(xs + w/2, R.sys_pct, width=w, color=C["B"], label="systemic")
for i, r in R.iterrows():
    ax2.text(xs[i], max(r.gut_pct, r.sys_pct) + 1.2,
             f"p={r.p:.1e}", ha="center", fontsize=6.5)
ax2.set_xticks(xs); ax2.set_xticklabels(R.subset, fontsize=7)
ax2.set_ylabel("TLR10$^+$ cells (%)")
ax2.set_title("b  Same difference inside each subset\n(not composition)", loc="left")
ax2.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.12))

ax3 = fig.add_subplot(gs[1, 1])
im = ax3.imshow(np.log10(hb.values + 1), cmap="Purples", aspect="auto")
ax3.set_xticks(range(len(GB))); ax3.set_xticklabels(GB, rotation=45, ha="right")
ax3.set_yticks(range(len(hb))); ax3.set_yticklabels(hb.index, fontsize=6.5)
ax3.set_title("c  The surface module is absent from\nepithelium; CASP4 is not (HPA, 154 types)", loc="left")
cb = fig.colorbar(im, ax=ax3, fraction=.045, pad=.03); cb.set_label("log$_{10}$(nCPM+1)", fontsize=7)
save(fig, f"{OUT}/Fig4_organ_placement")
