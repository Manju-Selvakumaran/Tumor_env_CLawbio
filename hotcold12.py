#!/usr/bin/env python
"""
Hot/Cold across 4 samples x 3 cancers on a COMMON (cross-cancer) scale.

- Pull absolute RSEM (batch-normalized) expression from cBioPortal for all
  panel genes and all 12 samples  -> comparable across cohorts.
- log2(x+1), then pooled z across the 12 tumors per gene -> one common scale.
- Layer 1 (infiltration, 10 genes) + Layer 2 (immunosuppression, 8 genes).
- Classify each tumor and summarise per-cancer spread + rank cancers.
"""
import json, urllib.request, statistics as st, os

base = "https://www.cbioportal.org/api"
HERE = r"C:\Users\manju\clawbio_hackathon"

LAYER1 = ["CD8A","CD4","CD3E","FOXP3","CD68","CD163","NCAM1","CD19","CD274","PDCD1"]
LAYER2 = ["TIGIT","LAG3","HAVCR2","PDCD1","CD274","IL10","TGFB1","CTLA4"]
GENES  = sorted(set(LAYER1)|set(LAYER2))
HI, LO = 0.5, -0.5

def post(path, body):
    req=urllib.request.Request(base+path,data=json.dumps(body).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json"})
    return json.load(urllib.request.urlopen(req,timeout=60))

samples = json.load(open(os.path.join(HERE,"hotcold_samples12.json")))

# map genes -> entrez
g = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", GENES)
sym2ez = {x["hugoGeneSymbol"]:x["entrezGeneId"] for x in g}
ez2sym = {v:k for k,v in sym2ez.items()}

import math
# fetch absolute RSEM per study, log2 transform
expr = {}   # (label, sample) -> {gene: log2 rsem}
for study, samps in samples.items():
    label = study.split("_")[0].upper()
    recs = post(f"/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data/fetch?projection=SUMMARY",
                {"sampleIds": samps, "entrezGeneIds": list(sym2ez.values())})
    tmp = {s: {} for s in samps}
    for r in recs:
        s = r["sampleId"]; gene = ez2sym.get(r["entrezGeneId"]); v = r.get("value")
        if gene and v is not None and s in tmp:
            tmp[s][gene] = math.log2(v + 1.0)
    for s in samps:
        expr[(label, s)] = tmp[s]

# pooled z across all 12 tumors, per gene (common scale)
keys = list(expr.keys())
gene_mean, gene_sd = {}, {}
for gene in GENES:
    vals = [expr[k][gene] for k in keys if gene in expr[k]]
    gene_mean[gene] = st.mean(vals)
    gene_sd[gene]   = st.pstdev(vals) or 1e-9
z = {k: {gene: (expr[k][gene]-gene_mean[gene])/gene_sd[gene]
         for gene in GENES if gene in expr[k]} for k in keys}

def mean(xs):
    xs=[x for x in xs if x is not None]; return sum(xs)/len(xs) if xs else None
def band(x): return "HIGH" if x>=HI else ("LOW" if x<=LO else "MID")
def classify(l1,l2):
    if l1>=HI: return "HOT - inflamed" if l2<HI else "HOT - inflamed but immunosuppressed"
    if l1<=LO: return "COLD - immune desert"
    return "IMMUNE-EXCLUDED / altered (suppressed)" if l2>=HI else "IMMUNE-EXCLUDED / altered"

rows=[]
for (label,s) in keys:
    l1=mean([z[(label,s)].get(gm) for gm in LAYER1])
    l2=mean([z[(label,s)].get(gm) for gm in LAYER2])
    rows.append({"cancer":label,"sample":s,"L1":round(l1,3),"L2":round(l2,3),
                 "band1":band(l1),"band2":band(l2),"call":classify(l1,l2),
                 "z":{g:round(z[(label,s)].get(g),2) for g in GENES if g in z[(label,s)]}})

# per-cancer summary (spread + cross-cancer rank on common scale)
summary={}
for lab in ["SKCM","CESC","PAAD"]:
    l1s=[r["L1"] for r in rows if r["cancer"]==lab]
    l2s=[r["L2"] for r in rows if r["cancer"]==lab]
    summary[lab]={"n":len(l1s),
                  "L1_mean":round(st.mean(l1s),3),"L1_min":min(l1s),"L1_max":max(l1s),
                  "L1_sd":round(st.pstdev(l1s),3),
                  "L2_mean":round(st.mean(l2s),3),"L2_min":min(l2s),"L2_max":max(l2s),
                  "L2_sd":round(st.pstdev(l2s),3)}

json.dump({"scale":"pooled z across 12 tumors (common, cross-cancer) from log2 batch-normalized RSEM",
           "layer1_panel":LAYER1,"layer2_panel":LAYER2,"thresholds":{"high":HI,"low":LO},
           "tumors":rows,"per_cancer":summary},
          open(os.path.join(HERE,"hotcold12_result.json"),"w"),indent=2)

# ---- print ----
print("Per-tumor (common cross-cancer scale):")
print(f"  {'Cancer':6}{'Sample':18}{'L1':>7}{'L2':>7}  Classification")
for lab in ["SKCM","CESC","PAAD"]:
    for r in [x for x in rows if x["cancer"]==lab]:
        print(f"  {r['cancer']:6}{r['sample']:18}{r['L1']:>7}{r['L2']:>7}  {r['call']}")
print("\nPer-cancer spread & cross-cancer comparison (Layer 1 infiltration):")
print(f"  {'Cancer':6}{'mean':>7}{'min':>7}{'max':>7}{'sd':>7}")
for lab,srt in sorted(summary.items(),key=lambda kv:-kv[1]['L1_mean']):
    print(f"  {lab:6}{srt['L1_mean']:>7}{srt['L1_min']:>7}{srt['L1_max']:>7}{srt['L1_sd']:>7}")
print("\nPer-cancer (Layer 2 immunosuppression):")
print(f"  {'Cancer':6}{'mean':>7}{'min':>7}{'max':>7}{'sd':>7}")
for lab,srt in sorted(summary.items(),key=lambda kv:-kv[1]['L2_mean']):
    print(f"  {lab:6}{srt['L2_mean']:>7}{srt['L2_min']:>7}{srt['L2_max']:>7}{srt['L2_sd']:>7}")
