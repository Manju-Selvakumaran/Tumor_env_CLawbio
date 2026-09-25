#!/usr/bin/env python
"""Render hot/cold visualizations from hotcold12_result.json (45 tumors)."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = r"C:\Users\manju\clawbio_hackathon"
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
res = json.load(open(os.path.join(HERE, "hotcold12_result.json")))
tumors = res["tumors"]

# palette (dataviz reference: first 3 categorical slots, validated all-pairs)
CANCER_COLOR = {"SKCM": "#2a78d6", "CESC": "#eb6834", "PAAD": "#1baf7a"}
CANCER_MARK  = {"SKCM": "o",       "CESC": "s",       "PAAD": "^"}       # relief: shape too
INK, SEC, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SURF = "#fcfcfb"
HI, LO = 0.5, -0.5

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "figure.facecolor": SURF, "axes.facecolor": SURF,
    "axes.edgecolor": BASE, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
})

# ----------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(figsize=(8.4, 7.0))
# threshold guides
for x in (LO, HI): ax.axvline(x, color=BASE, lw=1, ls="--", zorder=1)
ax.axhline(HI, color=BASE, lw=1, ls="--", zorder=1)
# region tints
ax.axvspan(-3, LO, color="#2a78d6", alpha=0.05, zorder=0)   # cold band
ax.axvspan(HI, 3,  color="#e34948", alpha=0.05, zorder=0)   # hot band
# region labels
ax.text(-2.35, 1.85, "COLD\nimmune desert", color=SEC, fontsize=10, ha="left", va="top", weight="bold")
ax.text(1.65, 1.85, "HOT\ninflamed", color=SEC, fontsize=10, ha="left", va="top", weight="bold")
ax.text(0.0, -2.3, "IMMUNE-EXCLUDED / altered", color=MUTED, fontsize=9, ha="center", va="bottom")
ax.text(2.55, 1.30, "+ suppressed\n(L2 high)", color=MUTED, fontsize=8, ha="right", va="center", style="italic")

for t in tumors:
    c = CANCER_COLOR[t["cancer"]]
    ax.scatter(t["L1"], t["L2"], s=70, c=c, marker=CANCER_MARK[t["cancer"]],
               edgecolors="white", linewidths=0.8, alpha=0.9, zorder=3)

ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5)
ax.set_xlabel("Layer 1 — immune infiltration (pooled z)", fontsize=11)
ax.set_ylabel("Layer 2 — immunosuppression (pooled z)", fontsize=11)
ax.set_title("Tumor immune phenotype — 45 TCGA tumors (common cross-cancer scale)",
             fontsize=12.5, weight="bold", pad=12)
ax.grid(True, color=GRID, lw=0.6, zorder=0)
for s in ("top", "right"): ax.spines[s].set_visible(False)
leg = [Line2D([0],[0], marker=CANCER_MARK[k], color="none", markerfacecolor=v,
              markeredgecolor="white", markersize=9,
              label=f"{k} (n={sum(1 for t in tumors if t['cancer']==k)})")
       for k, v in CANCER_COLOR.items()]
ax.legend(handles=leg, loc="lower right", frameon=False, fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "hotcold_scatter.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------- Figure 2
GENE_ORDER = ["CD8A","CD4","CD3E","FOXP3","CD68","CD163","NCAM1","CD19","CD274","PDCD1",
              "TIGIT","LAG3","HAVCR2","IL10","TGFB1","CTLA4"]
# order tumors: by cancer, then by L1 ascending (cold -> hot)
order = ["SKCM","CESC","PAAD"]
cols = sorted(tumors, key=lambda t: (order.index(t["cancer"]), t["L1"]))
mat = [[c["z"].get(g, 0.0) for c in cols] for g in GENE_ORDER]

diverging = LinearSegmentedColormap.from_list("bwr_ds", ["#184f95","#6da7ec","#f0efec","#eb8f7d","#c0392b"])
figh, axh = plt.subplots(figsize=(13.5, 6.2))
im = axh.imshow(mat, cmap=diverging, vmin=-2.2, vmax=2.2, aspect="auto")
axh.set_yticks(range(len(GENE_ORDER))); axh.set_yticklabels(GENE_ORDER, fontsize=9)
axh.set_xticks(range(len(cols)))
axh.set_xticklabels([c["sample"].replace("TCGA-","").rsplit("-",1)[0] for c in cols],
                    rotation=90, fontsize=6, color=MUTED)
# layer divider (after 10 Layer-1 genes)
axh.axhline(9.5, color=INK, lw=1.4)
axh.text(-3.0, 4.5, "Layer 1\ninfiltration", rotation=90, va="center", ha="center", fontsize=9, color=SEC, weight="bold")
axh.text(-3.0, 12.5, "Layer 2\nsuppression", rotation=90, va="center", ha="center", fontsize=9, color=SEC, weight="bold")
# cancer group separators + labels along top
b = 0
for lab in order:
    n = sum(1 for c in cols if c["cancer"] == lab)
    if b > 0: axh.axvline(b-0.5, color=INK, lw=1.6)
    axh.text(b + n/2 - 0.5, -1.1, f"{lab} (n={n})", ha="center", va="bottom",
             fontsize=10, weight="bold", color=CANCER_COLOR[lab])
    b += n
axh.set_title("Immune-marker expression across 45 tumors  (pooled z, ordered cold→hot within cancer)",
              fontsize=12.5, weight="bold", pad=26)
cbar = figh.colorbar(im, ax=axh, fraction=0.018, pad=0.01)
cbar.set_label("pooled z-score", fontsize=9); cbar.ax.tick_params(labelsize=8)
axh.set_xlabel("tumor sample", fontsize=10, color=SEC)
figh.tight_layout()
figh.savefig(os.path.join(FIG, "hotcold_heatmap.png"), dpi=200, bbox_inches="tight")
plt.close(figh)

print("Wrote:")
print(" ", os.path.join(FIG, "hotcold_scatter.png"))
print(" ", os.path.join(FIG, "hotcold_heatmap.png"))
