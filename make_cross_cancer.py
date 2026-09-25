#!/usr/bin/env python
"""
Cross-cancer reference matrix: score each cancer's tumors against EACH cancer's
own distribution (SKCM / CESC / PAAD reference), for Layer 1 and Layer 2.

cell(i, j) = mean layer score of cancer i's tumors measured on cancer j's yardstick
  z_g = (log2 RSEM of tumor - mean_j[g]) / sd_j[g]
  diagonal (i==j) ~ 0 by construction; off-diagonal = i relative to j's baseline.
"""
import json, os, math, statistics as st, urllib.request
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = r"C:\Users\manju\clawbio_hackathon"
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
base = "https://www.cbioportal.org/api"
def post(path, body):
    req = urllib.request.Request(base+path, data=json.dumps(body).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json"})
    return json.load(urllib.request.urlopen(req, timeout=60))

LAYER1 = ["CD8A","CD4","CD3E","FOXP3","CD68","CD163","NCAM1","CD19","CD274","PDCD1"]
LAYER2 = ["TIGIT","LAG3","HAVCR2","PDCD1","CD274","IL10","TGFB1","CTLA4"]
GENES  = sorted(set(LAYER1) | set(LAYER2))
CANCERS = ["SKCM", "CESC", "PAAD"]

samples = json.load(open(os.path.join(HERE, "hotcold_samples12.json")))
g = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", GENES)
sym2ez = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in g}; ez2sym = {v:k for k,v in sym2ez.items()}

# fetch absolute log2 RSEM for all tumors, keyed by cancer
expr = {c: [] for c in CANCERS}   # cancer -> list of {gene: log2}
for study, samps in samples.items():
    lab = study.split("_")[0].upper()
    recs = post(f"/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data/fetch?projection=SUMMARY",
                {"sampleIds": samps, "entrezGeneIds": list(sym2ez.values())})
    per = {s: {} for s in samps}
    for r in recs:
        gene = ez2sym.get(r["entrezGeneId"]); v = r.get("value")
        if gene and v is not None: per[r["sampleId"]][gene] = math.log2(v + 1.0)
    expr[lab] = list(per.values())

# reference stats per cancer: gene -> (mean, sd)
ref = {}
for c in CANCERS:
    ref[c] = {}
    for gene in GENES:
        vals = [t[gene] for t in expr[c] if gene in t]
        ref[c][gene] = (st.mean(vals), st.pstdev(vals) or 1e-9)

def layer_score(tumor, refc, panel):
    zs = [(tumor[gm] - ref[refc][gm][0]) / ref[refc][gm][1] for gm in panel if gm in tumor]
    return sum(zs) / len(zs) if zs else 0.0

# 3x3 matrices: row = tumor's cancer i, col = reference cancer j
def matrix(panel):
    return [[st.mean([layer_score(t, j, panel) for t in expr[i]]) for j in CANCERS]
            for i in CANCERS]
M1 = matrix(LAYER1); M2 = matrix(LAYER2)

# --- plot ---
INK, SEC, MUTED, SURF = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Segoe UI","DejaVu Sans"],
                     "figure.facecolor":SURF,"axes.facecolor":SURF,"text.color":INK})
div = LinearSegmentedColormap.from_list("bwr", ["#184f95","#6da7ec","#f0efec","#eb8f7d","#c0392b"])
vlim = max(abs(v) for M in (M1, M2) for row in M for v in row)
vlim = max(vlim, 0.5)

fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2))
for ax, (M, name) in zip(axes, [(M1, "Layer 1 — infiltration"), (M2, "Layer 2 — immunosuppression")]):
    im = ax.imshow(M, cmap=div, vmin=-vlim, vmax=vlim, aspect="equal")
    ax.set_xticks(range(3)); ax.set_xticklabels([f"{c}\nreference" for c in CANCERS], fontsize=9)
    ax.set_yticks(range(3)); ax.set_yticklabels(CANCERS, fontsize=10)
    for i in range(3):
        for j in range(3):
            v = M[i][j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=11,
                    color="white" if abs(v) > 0.55*vlim else INK, weight="bold")
            if i == j:
                ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                             edgecolor=MUTED, lw=1.4, ls="--"))
    ax.set_title(name, fontsize=11.5, weight="bold", pad=8)
    ax.set_xlabel("scored against  (reference cancer's yardstick)", fontsize=9, color=SEC)
axes[0].set_ylabel("tumors from  (cancer type)", fontsize=9, color=SEC)
cbar = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.03)
cbar.set_label("mean z vs reference   (+ hotter / − colder than that reference)", fontsize=9)
cbar.ax.tick_params(labelsize=8)
fig.suptitle("Cross-cancer comparison — each cancer scored on every cancer's reference distribution",
             fontsize=13, weight="bold", y=1.02)
out = os.path.join(FIG, "hotcold_cross_cancer.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)

json.dump({"cancers": CANCERS, "layer1_matrix": M1, "layer2_matrix": M2,
           "note": "cell[i][j] = mean layer z of cancer i tumors on cancer j reference"},
          open(os.path.join(HERE, "hotcold_cross_cancer.json"), "w"), indent=2)

print("Layer 1 (rows=tumor cancer, cols=reference):")
print("        " + "  ".join(f"{c:>7}" for c in CANCERS))
for i, c in enumerate(CANCERS):
    print(f"{c:6} " + "  ".join(f"{M1[i][j]:>+7.2f}" for j in range(3)))
print("Wrote", out)
