"""共通スタイル＋入出力パス。図は png/pdf/svg の3形式で出す。seed 固定・生データ read-only。

パス規約:
  DATA_ROOT  入力データの置き場。環境変数 DATA_ROOT で上書きできる。
             既定は開発機の /Volumes/M4_SSD/projects。
             どのファイルをどこに置くかは README.md の Data manifest を見ること。
  FIGDIR     このファイルのあるディレクトリ。図と fig*_source.csv の出力先。
"""
import os, matplotlib as mpl, matplotlib.pyplot as plt, numpy as np

DATA_ROOT = os.environ.get("DATA_ROOT", "/Volumes/M4_SSD/projects")
FIGDIR = os.path.dirname(os.path.abspath(__file__))


def need(relpath):
    """DATA_ROOT 以下の入力パスを返す。無ければ、何をどこに置くかを言って止まる。"""
    p = os.path.join(DATA_ROOT, relpath)
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        raise SystemExit(
            f"入力が見つからない（または空）: {p}\n"
            f"  DATA_ROOT={DATA_ROOT}\n"
            f"  README.md の Data manifest を見て {relpath} を配置するか、"
            f"DATA_ROOT を設定してください。\n"
            f"  ※ 0 バイトの場合は外付けディスクへのシンボリックリンクが切れている可能性があります。"
        )
    return p


mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif", "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
})
# colour-blind safe
C = {"B": "#0072B2", "mye": "#D55E00", "nk": "#009E73", "epi": "#CC79A7",
     "plasma": "#8C6D31", "stroma": "#999999", "gran": "#E69F00",
     "up": "#B2182B", "down": "#2166AC", "flat": "#8A8A8A"}
SEED = 0
np.random.seed(SEED)

def save(fig, path_noext):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{path_noext}.{ext}")
    print("wrote", path_noext + ".{png,pdf,svg}")
