#!/usr/bin/env python3
"""
Fig 3: what moves is the gate, not the receiver.
(a) SCP3755 患者内ペア Δ（GBP2/STAT1/IRF1 vs CASP4 vs TLR10）
(b) GSE214695 患者レベル GBP2（HC vs IBD, 細胞種別）
(c) CASP4+GBP2 二重陽性率（licensing）の拡大
"""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import C, save, DATA_ROOT, FIGDIR
from scipy.stats import wilcoxon, mannwhitneyu

D = DATA_ROOT; OUT = FIGDIR
LINS = ["Myeloid", "B", "T_NK_ILC", "Epithelium", "Plasma"]
GEN = ["GBP2", "STAT1", "IRF1", "CASP4", "TLR10"]

# ---- (a) SCP3755 paired ----
E = pd.read_pickle(f"{D}/scp3755_uc_atlas/E_logcpm.pkl")
E = E[E.n_cells >= 20]
def paired(gene, lineage):
    d = E[E.lin == lineage]; out = []
    for p, gp in d.groupby("PatientID"):
        for ct, g2 in gp.groupby("celltype"):
            i = g2[g2.Biopsy_INFLAMMATION == "Inflamed"][gene]
            u = g2[g2.Biopsy_INFLAMMATION == "Uninflamed"][gene]
            if len(i) and len(u): out.append((p, i.mean() - u.mean()))
    return pd.DataFrame(out, columns=["pt", "d"]).groupby("pt").d.mean()

# 発現量マッチ帰無（07_gbp2_paired.py と同一手順、seed=0、シグネチャ遺伝子は除外）
PB = f"{D}/scp3755_uc_atlas/pb_counts.csv"
hdr = open(PB).readline().rstrip("\n").split(",")[1:]
names, tot = [], np.zeros(len(hdr))
for line in open(PB):
    g, _, rest = line.partition(",")
    if g == "" or g == hdr[0] and not names and rest.startswith("SAM"):
        continue
    try:
        v = np.array(rest.rstrip("\n").split(","), dtype=float)
    except ValueError:
        continue
    names.append(g); tot += v
rng = np.random.default_rng(0)
pool = set(rng.choice(np.array(names), 1200, replace=False).tolist())
EXCL = {g for g in pool if g.startswith("GBP") or g.startswith("CASP")} | set(GEN) | {"S100A8"}
pool = sorted(pool - EXCL)
vals = {}
need = set(pool)
for line in open(PB):
    g, _, rest = line.partition(",")
    if g in need:
        try:
            vals[g] = np.log1p(np.array(rest.rstrip("\n").split(","), dtype=float) / tot * 1e6)
        except ValueError:
            pass
        need.discard(g)
POOL = pd.DataFrame(vals, index=hdr).join(
    E[["PatientID", "celltype", "Biopsy_INFLAMMATION", "lin"]], how="inner")
print("pool genes:", POOL.shape)

def paired_pool(gene, lineage):
    d = POOL[POOL.lin == lineage]; out = []
    for ptid, gp in d.groupby("PatientID"):
        v = []
        for ct, g2 in gp.groupby("celltype"):
            i = g2[g2.Biopsy_INFLAMMATION == "Inflamed"][gene]
            u = g2[g2.Biopsy_INFLAMMATION == "Uninflamed"][gene]
            if len(i) and len(u): v.append(i.mean() - u.mean())
        if v: out.append(np.mean(v))
    return np.mean(out) if len(out) >= 5 else np.nan

M = pd.DataFrame(index=LINS, columns=GEN, dtype=float)
P = M.copy(); NP = M.copy()
for ln in LINS:
    base = E[E.lin == ln][GEN].mean()
    pbase = POOL[pool][POOL.lin == ln].mean()
    for g in GEN:
        v = paired(g, ln)
        M.loc[ln, g] = v.mean()
        P.loc[ln, g] = wilcoxon(v.values).pvalue if len(v) >= 6 else np.nan
        cand = [x for x in pbase.index if abs(pbase[x] - base[g]) < 0.15][:120]
        null = [abs(paired_pool(x, ln)) for x in cand]
        null = [x for x in null if not np.isnan(x)]
        NP.loc[ln, g] = float(np.mean(np.array(null) >= abs(M.loc[ln, g]))) if null else np.nan
M.to_csv(f"{OUT}/fig3a_delta.csv"); P.to_csv(f"{OUT}/fig3a_p.csv"); NP.to_csv(f"{OUT}/fig3a_nullp.csv")
print(M.round(3).to_string()); print("raw p\n" + P.round(4).to_string()); print("null-calibrated p\n" + NP.round(3).to_string())

# ---- (b,c) GSE214695 ----
pl = pd.read_csv(f"{D}/gbp2_cytosolic_arm/patient_lineage.csv")
pl["grp"] = np.where(pl.group == "HC", "HC", "IBD")
KEEP = ["Myeloid", "B", "Epithelium", "T_NK_ILC", "Stromal", "Plasma"]
pl = pl[pl.lin.isin(KEEP)]

# ---- draw ----
fig = plt.figure(figsize=(7.2, 5.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.25, 1.0],
                      hspace=0.55, wspace=0.35)

# (a) dot matrix
ax = fig.add_subplot(gs[0, :])
for i, g in enumerate(GEN):
    for j, ln in enumerate(LINS):
        d = M.loc[ln, g]; p = P.loc[ln, g]
        col = C["up"] if d > 0 else C["down"]
        npv = NP.loc[ln, g]
        passed = (p < 0.05) and (npv < 0.05)
        ax.scatter(j, -i, s=40 + 260 * min(abs(d), 1.5) / 1.5,
                   facecolor=(col if passed else "none"),
                   edgecolor=col, linewidth=(0 if passed else 1.1), alpha=.85)
        if passed:
            ax.text(j, -i - 0.30, "*", ha="center", va="center", fontsize=11, color="k")
ax.set_xticks(range(len(LINS))); ax.set_xticklabels([l.replace("_", "/") for l in LINS])
ax.set_yticks([-i for i in range(len(GEN))]); ax.set_yticklabels(GEN)
ax.set_xlim(-.6, len(LINS) - .4); ax.set_ylim(-len(GEN) + .4, .7)
ax.set_title("a  Within-patient inflamed $-$ uninflamed (SCP3755). "
             "Red = up, blue = down; size = |$\\Delta$|; filled + * = passes expression-matched null", loc="left")
for lab, xx in (("IFN-inducible gate", -0.45), ):
    pass
ax.annotate("gate", xy=(len(LINS) - .45, -1), xytext=(len(LINS) - .1, -1),
            fontsize=7.5, va="center", color="#444")
ax.annotate("receiver", xy=(len(LINS) - .45, -3), xytext=(len(LINS) - .1, -3),
            fontsize=7.5, va="center", color="#444")
ax.annotate("reader", xy=(len(LINS) - .45, -4), xytext=(len(LINS) - .1, -4),
            fontsize=7.5, va="center", color="#444")

# (b) GBP2 patient-level
ax2 = fig.add_subplot(gs[1, 0])
xs = np.arange(len(KEEP))
for k, grp in enumerate(["HC", "IBD"]):
    for j, ln in enumerate(KEEP):
        v = pl[(pl.lin == ln) & (pl.grp == grp)].GBP2.values
        x = xs[j] + (k - .5) * .34
        ax2.scatter(np.full(len(v), x) + np.random.uniform(-.045, .045, len(v)), v,
                    s=10, color=(C["flat"] if grp == "HC" else C["mye"]), alpha=.7, lw=0)
        ax2.plot([x - .12, x + .12], [np.mean(v)] * 2, lw=2,
                 color=("#444" if grp == "HC" else C["mye"]))
ax2.set_xticks(xs); ax2.set_xticklabels([l.replace("_", "/") for l in KEEP], rotation=30, ha="right")
ax2.set_ylabel("GBP2 (patient mean, log-norm)")
ax2.set_title("b  GBP2 rises lineage-selectively\n(GSE214695: HC n=6, IBD n=12)", loc="left")
ax2.scatter([], [], s=14, color=C["flat"], label="HC"); ax2.scatter([], [], s=14, color=C["mye"], label="IBD")
ax2.legend(loc="upper left", ncol=2)

# (c) licensing
ax3 = fig.add_subplot(gs[1, 1])
sub = ["Myeloid", "B", "Stromal", "Plasma"]
w = .34
for k, grp in enumerate(["HC", "IBD"]):
    vals = [pl[(pl.lin == ln) & (pl.grp == grp)].co_CASP4_GBP2.mean() for ln in sub]
    ax3.bar(np.arange(len(sub)) + (k - .5) * w, vals, width=w,
            color=(C["flat"] if grp == "HC" else C["mye"]), label=grp)
ax3.set_xticks(range(len(sub))); ax3.set_xticklabels(sub, rotation=30, ha="right")
ax3.set_ylabel("CASP4$^+$GBP2$^+$ cells (%)")
ax3.set_title("c  Licensed fraction expands\nin myeloid and B, not stroma/plasma", loc="left")
ax3.legend(loc="upper right", ncol=2)
save(fig, f"{OUT}/Fig3_gate_not_receiver")

# 参考: 患者レベル検定（本文用）
for ln in KEEP:
    a = pl[(pl.lin == ln) & (pl.grp == "HC")]; b = pl[(pl.lin == ln) & (pl.grp == "IBD")]
    p1 = mannwhitneyu(a.GBP2, b.GBP2).pvalue
    p2 = mannwhitneyu(a.co_CASP4_GBP2, b.co_CASP4_GBP2).pvalue
    print(f"{ln:11s} GBP2 {a.GBP2.mean():.3f}->{b.GBP2.mean():.3f} p={p1:.4f} | "
          f"co {a.co_CASP4_GBP2.mean():.1f}%->{b.co_CASP4_GBP2.mean():.1f}% p={p2:.4f}")
