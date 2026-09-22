#!/usr/bin/env python3
"""
Fig 2: TLR10 is a fixed installation.
(a) 活性化で動かない: GSE327707 分取B細胞 6ドナー paired 48h, 4刺激。内部対照 LY86/TLR1。
(b) 炎症で動かない: 3コホート2プラットフォームの TLR10 Δ。
生データ read-only / seed 固定。
"""
import os, sys, gzip, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import C, save, DATA_ROOT, FIGDIR
from scipy.stats import wilcoxon

D = DATA_ROOT
OUT = FIGDIR

# ---------- (a) GSE327707 ----------
G = f"{D}/bcell_lipidA_whitespace/GSE327707"
cd = pd.read_csv(f"{G}/data/coldata.csv")
X = pd.read_csv(f"{G}/data/GSE327707_raw_counts.csv.gz", index_col=0)
if not X.index.astype(str).str.startswith("ENSG").any():
    sym = X.index.astype(str)
else:
    ref = pd.read_csv(f"{D}/spencer_galt/geo/GSE327417/extracted/GSM9656505_318B-2_features.tsv.gz",
                      sep="\t", header=None, names=["ens", "symbol", "type"])
    m = dict(zip(ref.ens, ref.symbol))
    sym = pd.Index([m.get(str(i).split(".")[0], str(i)) for i in X.index])
X.index = sym
cpm = X / X.sum() * 1e6
L = np.log2(cpm + 1)

GENES = ["TLR10", "LY86", "TLR1"]
STIMS = ["LPS", "IFNa", "Cytostim", "HEK BCS (CD79B)"]
LAB = {"LPS": "LPS", "IFNa": "IFNα", "Cytostim": "Cytostim", "HEK BCS (CD79B)": "BCR"}
b = cd[cd.celltype == "B_cell"]
rows = []
for g in GENES:
    for s in STIMS:
        for dn in b.donor.unique():
            t = b[(b.donor == dn) & (b.treatment == s)].lib
            u = b[(b.donor == dn) & (b.treatment == "Baseline")].lib
            t = [c for c in t if c in L.columns]; u = [c for c in u if c in L.columns]
            if t and u:
                rows.append((g, s, dn, L.loc[g, t].mean() - L.loc[g, u].mean()))
A = pd.DataFrame(rows, columns=["gene", "stim", "donor", "d"])
stats = (A.groupby(["gene", "stim"]).d
           .agg(["mean", "count"]).reset_index())
pv = []
for _, r in stats.iterrows():
    v = A[(A.gene == r.gene) & (A.stim == r.stim)].d.values
    pv.append(wilcoxon(v).pvalue if len(v) >= 6 else np.nan)
stats["p"] = pv
stats.to_csv(f"{OUT}/fig2a_source.csv", index=False)
print(stats.to_string())

# ---------- (b) three cohorts ----------
# GSE214695 (10x, cross-sectional HC vs IBD) : patient-level means
pl = pd.read_csv(f"{D}/gbp2_cytosolic_arm/patient_lineage.csv")
pl["grp"] = np.where(pl.group == "HC", "HC", "IBD")
g1 = (pl.groupby(["lin", "grp"]).TLR10.mean().unstack())
g1["d"] = g1["IBD"] - g1["HC"]

# SCP3755 (BD Rhapsody, within-patient pairs)
E = pd.read_pickle(f"{D}/scp3755_uc_atlas/E_logcpm.pkl")
pm = pd.read_csv(f"{D}/scp3755_uc_atlas/pb_meta.csv").set_index("group")
sm = (pd.read_csv(f"{D}/scp3755_uc_atlas/scp3755_metadata.csv",
                  usecols=["SAMID", "PatientID", "Biopsy_INFLAMMATION"])
        .drop_duplicates("SAMID").set_index("SAMID"))
for c, src in (("SAMID", pm[["SAMID", "celltype", "n_cells"]]), ("PatientID", None)):
    pass
if "SAMID" not in E.columns:
    E = E.join(pm[["SAMID", "celltype", "n_cells"]])
if "PatientID" not in E.columns:
    E = E.join(sm, on="SAMID")

def lin(c):
    if c in ("Macrophage", "Resident macrophage", "Monocyte-derived macrophage", "Dendritic cell"): return "Myeloid"
    if c.startswith("Neutrophil") or c == "Eosinophil": return "Granulocyte"
    if c in ("Follicular B", "GC", "Cycling B"): return "B"
    if c == "Plasma": return "Plasma"
    if c in ("CD4 T Cell", "CD8 T Cell", "Cycling T Cell", "NK", "ILC"): return "T_NK_ILC"
    if c in ("Enterocytes", "Enterocytes BEST4", "Goblet", "Stem Cells", "TA", "Tuft", "Enteroendocrine"): return "Epithelium"
    if "fibroblast" in c or c in ("Pericyte", "SMC", "Glial"): return "Stromal"
    if "endothelial" in c or c == "Lymphatics": return "Endothelium"
    return "Other"
E["lin"] = E.celltype.map(lin)
E = E[E.n_cells >= 20]

def paired(gene, lineage):
    d = E[E.lin == lineage]; out = []
    for p, gp in d.groupby("PatientID"):
        for ct, g2 in gp.groupby("celltype"):
            i = g2[g2.Biopsy_INFLAMMATION == "Inflamed"][gene]
            u = g2[g2.Biopsy_INFLAMMATION == "Uninflamed"][gene]
            if len(i) >= 1 and len(u) >= 1:
                out.append((p, ct, i.mean() - u.mean()))
    return pd.DataFrame(out, columns=["pt", "ct", "d"])

LINS = ["Myeloid", "B", "T_NK_ILC", "Epithelium", "Plasma"]
rows = []
for ln in LINS:
    for gene in ["TLR10", "S100A8"]:
        dd = paired(gene, ln)
        if len(dd) == 0: continue
        pt = dd.groupby("pt").d.mean()
        p = wilcoxon(pt.values).pvalue if len(pt) >= 6 else np.nan
        rows.append((ln, gene, pt.mean(), len(pt), int((dd.d > 0).sum()), len(dd), p))
B3 = pd.DataFrame(rows, columns=["lin", "gene", "d", "n_pt", "up", "n_unit", "p"])
B3.to_csv(f"{OUT}/fig2b_source.csv", index=False)
print(B3.to_string())

# ---------- draw ----------
fig = plt.figure(figsize=(7.2, 3.1))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.42)

ax = fig.add_subplot(gs[0])
xs = np.arange(len(STIMS)); w = 0.26
col = {"TLR10": C["B"], "LY86": C["mye"], "TLR1": C["nk"]}
for k, g in enumerate(GENES):
    sub = A[A.gene == g]
    for j, s in enumerate(STIMS):
        v = sub[sub.stim == s].d.values
        ax.scatter(np.full(len(v), xs[j] + (k - 1) * w) + np.random.uniform(-.05, .05, len(v)),
                   v, s=9, color=col[g], alpha=.55, lw=0)
        ax.plot([xs[j] + (k - 1) * w - .1, xs[j] + (k - 1) * w + .1], [v.mean()] * 2,
                color=col[g], lw=2)
ax.axhline(0, color="k", lw=.6)
ax.set_xticks(xs); ax.set_xticklabels([LAB[s] for s in STIMS])
ax.set_ylabel("paired $\\Delta$ log$_2$ CPM (48 h $-$ baseline)")
ax.set_title("a  Sorted human B cells do not move TLR10\n(GSE327707, 6 donors, paired)", loc="left")
for g in GENES:
    ax.scatter([], [], s=14, color=col[g], label=g)
ax.legend(loc="lower left", ncol=3, handletextpad=.2, columnspacing=.9)

ax2 = fig.add_subplot(gs[1])
t = B3[B3.gene == "TLR10"].set_index("lin").reindex(LINS)
s8 = B3[B3.gene == "S100A8"].set_index("lin").reindex(LINS)
y = np.arange(len(LINS))[::-1]
ax2.barh(y + .18, s8.d.values, height=.34, color=C["flat"], alpha=.5, label="S100A8 (positive control)")
ax2.barh(y - .18, t.d.values, height=.34, color=C["B"], label="TLR10")
for i, ln in enumerate(LINS):
    r = t.loc[ln]
    ax2.text(0.06 if r.d < 0 else r.d + .06, y[i] - .18,
             f"{r.up}/{r.n_unit} units up", va="center", fontsize=6.5, color="#333")
ax2.axvline(0, color="k", lw=.6)
ax2.set_yticks(y); ax2.set_yticklabels([l.replace("_", "/") for l in LINS])
ax2.set_xlabel("within-patient $\\Delta$ log CPM (inflamed $-$ uninflamed)")
ax2.set_title("b  TLR10 does not move with inflammation\n(SCP3755, BD Rhapsody, 14–17 patients)", loc="left")
ax2.set_xlim(t.d.min()-0.45, max(s8.d.max(), 0)+0.55)
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
save(fig, f"{OUT}/Fig2_fixed_installation")
