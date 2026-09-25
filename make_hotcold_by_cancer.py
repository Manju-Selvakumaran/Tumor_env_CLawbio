#!/usr/bin/env python
"""Cross-cancer hot/cold heatmap: cancer type x phenotype category (fraction of tumors)."""
import json, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = r"C:\Users\manju\clawbio_hackathon"
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
tumors = json.load(open(os.path.join(HERE, "hotcold12_result.json")))["tumors"]

CATS = ["Cold", "Excluded", "Excluded\n(suppressed)", "Hot\n(inflamed)", "Hot\n(suppressed)"]
def cat(call):
    if "COLD" in call: return "Cold"
    if "HOT - inflamed but" in call: return "Hot\n(suppressed)"
    if "HOT - inflamed" in call: return "Hot\n(inflamed)"
    if "suppressed" in call: return "Excluded\n(suppressed)"
    return "Excluded"

rows = ["SKCM", "CESC", "PAAD"]
counts = {r: {c: 0 for c in CATS} for r in rows}
for t in tumors:
    counts[t["cancer"]][cat(t.get("classification") or t.get("call"))] += 1

frac = [[counts[r][c] / sum(counts[r].values()) for c in CATS] for r in rows]

INK, SEC, MUTED, SURF = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
CANCER_COLOR = {"SKCM": "#2a78d6", "CESC": "#eb6834", "PAAD": "#1baf7a"}
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK})
# cold->hot diverging ramp for the column axis feel; cell shade = fraction (sequential)
seq = LinearSegmentedColormap.from_list("frac", ["#f4f8ec", "#cfe6c6", "#7bc47f",
                                                 "#eda100", "#eb6834", "#c0392b"])

fig, ax = plt.subplots(figsize=(9.2, 4.4))
im = ax.imshow(frac, cmap=seq, vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(CATS))); ax.set_xticklabels(CATS, fontsize=10)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([f"{r} (n={sum(counts[r].values())})" for r in rows], fontsize=11)
for yt, r in zip(ax.get_yticklabels(), rows): yt.set_color(CANCER_COLOR[r])

# annotate each cell with count + percentage
for i, r in enumerate(rows):
    tot = sum(counts[r].values())
    for j, c in enumerate(CATS):
        n = counts[r][c]; pct = 100 * n / tot
        txt = f"{n}\n{pct:.0f}%"
        ax.text(j, i, txt, ha="center", va="center", fontsize=10,
                color="white" if frac[i][j] > 0.45 else INK, weight="bold")

# arrow annotation cold -> hot
ax.annotate("", xy=(4.5, -0.85), xytext=(-0.5, -0.85), annotation_clip=False,
            arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.4))
ax.text(-0.5, -1.05, "COLD", color="#2a78d6", fontsize=9, weight="bold", ha="left", va="bottom", clip_on=False)
ax.text(4.5, -1.05, "HOT", color="#c0392b", fontsize=9, weight="bold", ha="right", va="bottom", clip_on=False)

ax.set_title("Hot vs Cold immune phenotype across cancer types  (fraction of tumors, common scale)",
             fontsize=12.5, weight="bold", pad=30)
ax.set_xticks([x - 0.5 for x in range(1, len(CATS))], minor=True)
ax.set_yticks([y - 0.5 for y in range(1, len(rows))], minor=True)
ax.grid(which="minor", color=SURF, lw=3)
ax.tick_params(which="minor", length=0)
cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("fraction of the cancer's tumors", fontsize=9); cbar.ax.tick_params(labelsize=8)
fig.tight_layout()
out = os.path.join(FIG, "hotcold_by_cancer_heatmap.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)

# console
print(f"{'':6} " + " ".join(f"{c.replace(chr(10),' '):>22}" for c in CATS))
for r in rows:
    print(f"{r:6} " + " ".join(f"{counts[r][c]:>22}" for c in CATS))
print("Wrote", out)
