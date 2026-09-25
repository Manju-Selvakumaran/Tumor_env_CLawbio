#!/usr/bin/env python
"""Expression heatmap (absolute log2 RSEM) for the 45 tumors, grouped by cancer."""
import json, os, math, urllib.request
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

GENE_ORDER = ["CD8A","CD4","CD3E","FOXP3","CD68","CD163","NCAM1","CD19","CD274","PDCD1",
              "TIGIT","LAG3","HAVCR2","IL10","TGFB1","CTLA4"]
samples = json.load(open(os.path.join(HERE,"hotcold_samples12.json")))
res = {t["sample"]: t for t in json.load(open(os.path.join(HERE,"hotcold12_result.json")))["tumors"]}

g = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", GENE_ORDER)
sym2ez = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in g}; ez2sym = {v:k for k,v in sym2ez.items()}

# fetch absolute RSEM -> log2
log2 = {}   # sample -> {gene: log2}
cancer_of = {}
for study, samps in samples.items():
    lab = study.split("_")[0].upper()
    recs = post(f"/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data/fetch?projection=SUMMARY",
                {"sampleIds": samps, "entrezGeneIds": list(sym2ez.values())})
    per = {s: {} for s in samps}
    for r in recs:
        gene = ez2sym.get(r["entrezGeneId"]); v = r.get("value")
        if gene and v is not None: per[r["sampleId"]][gene] = math.log2(v + 1.0)
    for s in samps:
        log2[s] = per[s]; cancer_of[s] = lab

# order: cancer, then by Layer-1 z (cold -> hot) using existing results
order = ["SKCM","CESC","PAAD"]
cols = sorted(log2, key=lambda s: (order.index(cancer_of[s]), res[s]["L1"]))
mat = [[log2[s].get(gene, 0.0) for s in cols] for gene in GENE_ORDER]

INK, SEC, MUTED, BASE, SURF = "#0b0b0b", "#52514e", "#898781", "#c3c2b7", "#fcfcfb"
CANCER_COLOR = {"SKCM":"#2a78d6","CESC":"#eb6834","PAAD":"#1baf7a"}
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Segoe UI","DejaVu Sans"],
                     "figure.facecolor":SURF,"axes.facecolor":SURF,"text.color":INK})
# sequential single-hue blue ramp (dataviz reference: light->dark = low->high)
seq = LinearSegmentedColormap.from_list("blue_seq",
        ["#eef5fd","#cde2fb","#86b6ef","#3987e5","#256abf","#0d366b"])

fig, ax = plt.subplots(figsize=(13.5, 6.2))
vmax = max(max(r) for r in mat)
im = ax.imshow(mat, cmap=seq, vmin=0, vmax=vmax, aspect="auto")
ax.set_yticks(range(len(GENE_ORDER))); ax.set_yticklabels(GENE_ORDER, fontsize=9)
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([s.replace("TCGA-","").rsplit("-",1)[0] for s in cols],
                   rotation=90, fontsize=6, color=MUTED)
ax.axhline(9.5, color=INK, lw=1.4)
ax.text(-3.0, 4.5, "Layer 1\ninfiltration", rotation=90, va="center", ha="center",
        fontsize=9, color=SEC, weight="bold")
ax.text(-3.0, 12.5, "Layer 2\nsuppression", rotation=90, va="center", ha="center",
        fontsize=9, color=SEC, weight="bold")
b = 0
for lab in order:
    n = sum(1 for s in cols if cancer_of[s] == lab)
    if b > 0: ax.axvline(b-0.5, color=INK, lw=1.6)
    ax.text(b + n/2 - 0.5, -1.1, f"{lab} (n={n})", ha="center", va="bottom",
            fontsize=10, weight="bold", color=CANCER_COLOR[lab])
    b += n
ax.set_title("Immune-marker EXPRESSION across 45 tumors  (log2 RSEM+1, ordered cold\u2192hot within cancer)",
             fontsize=12.5, weight="bold", pad=26)
cbar = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01)
cbar.set_label("log2(RSEM+1) expression", fontsize=9); cbar.ax.tick_params(labelsize=8)
ax.set_xlabel("tumor sample", fontsize=10, color=SEC)
fig.tight_layout()
out = os.path.join(FIG, "hotcold_expression_heatmap.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
print("Wrote", out)
