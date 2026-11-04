#!/usr/bin/env python3
"""Fig 1c dataset 1 (GSE266546, Crohn's terminal ileum whole tissue, 110,343 cells).
生データから系統別平均→遺伝子ごとに系統間 z 化を再計算し、発現マッチ帰無で較正する。
2026-09-27 監査で「この列の計算元がリポジトリに無い」と判明したため復元。"""
import gzip, os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA_ROOT
D = os.path.join(DATA_ROOT, "GSE266546")
LMAP={'B':'B','B Proliferating':'B',
      'Mast':'Myeloid','Neutrophils':'Myeloid','Recruited Macro.':'Myeloid','Resident Macro.':'Myeloid',
      'NK':'NK','T':'T','T Proliferating':'T',
      'Plasma':'Plasma','Plasma Proliferating':'Plasma'}
LN=["B","Myeloid","NK","T","Plasma"]
ann=pd.concat([pd.read_csv(f"{D}/GSE266546_Annotation_TI_{k}.txt.gz",sep="\t") for k in ("Endoscopy","Surgical")])
lin=dict(zip(ann.Barcodes.astype(str),ann.Celltypes.map(LMAP)))
print(f"cells annotated {len(lin):,}")
acc=None; ncell=np.zeros(len(LN))
for k in ("Endoscopy","Surgical"):
    f=f"{D}/GSE266546_Raw_counts_TI_{k}.txt.gz"
    cols=pd.read_csv(f,sep="\t",nrows=0).columns.tolist()
    lab=np.array([LN.index(lin.get(c)) if lin.get(c) in LN else -1 for c in cols])
    # pass 1: per-cell totals
    tot=np.zeros(len(cols))
    for ch in pd.read_csv(f,sep="\t",chunksize=200,index_col=0):
        tot+=ch.values.sum(axis=0)
    tot[tot==0]=1
    for i,l in enumerate(LN): ncell[i]+=int((lab==i).sum())
    print(f"  {k}: {len(cols):,} cells; per-lineage {[int((lab==i).sum()) for i in range(len(LN))]}")
    # pass 2: log1p(CPM) sums per lineage
    parts=[]
    for ch in pd.read_csv(f,sep="\t",chunksize=200,index_col=0):
        v=np.log1p(ch.values/tot*1e6)
        parts.append(pd.DataFrame({l:v[:,lab==i].sum(axis=1) for i,l in enumerate(LN)},index=ch.index))
    s=pd.concat(parts)
    acc=s if acc is None else acc.add(s,fill_value=0)
M=acc.div(ncell,axis=1)
M=M[~M.index.duplicated(keep='first')]
M.to_csv("fig1c_dataset1_lineage_means.csv")
z=M.sub(M.mean(axis=1),axis=0).div(M.std(axis=1,ddof=0)+1e-9,axis=0)
overall=M.mean(axis=1); expressed=overall[overall>0.005].index
print(f"\ngenes {len(M):,}  expressed {len(expressed):,}  cells/lineage {dict(zip(LN,ncell.astype(int)))}")
MS={'TLR10':2.5,'CD180':1.9,'LY86':2.3,'AOAH':None,'TLR4':None,'LY96':None}
print(f"{'gene':7s} {'B z':>7s} {'manuscript':>11s} {'null n':>8s} {'%ile':>8s}  verdict")
rows=[]
for g in ['TLR10','CD180','LY86','AOAH','TLR4','LY96']:
    if g not in M.index: print(f"{g:7s} ABSENT"); continue
    o=overall[g]; pool=[x for x in expressed if 0.9*o<=overall[x]<=1.1*o and x!=g]
    pct=(z.loc[pool,'B']<z.loc[g,'B']).mean()
    v='** specific' if pct>0.95 else 'INSIDE NULL'
    print(f"{g:7s} {z.loc[g,'B']:7.2f} {str(MS[g] or '-'):>11s} {len(pool):8d} {pct:7.1%}  {v}")
    rows.append(dict(gene=g,B_z=z.loc[g,'B'],null_n=len(pool),pctile=pct,verdict=v))
pd.DataFrame(rows).to_csv("fig1c_dataset1_null.csv",index=False)
z.loc[['TLR10','CD180','LY86','AOAH','TLR4','LY96']].to_csv("fig1c_dataset1_z.csv")
print("\nfull z table:\n"+z.loc[[g for g in ['TLR10','CD180','LY86','AOAH','TLR4','LY96'] if g in z.index]].round(2).to_string())
