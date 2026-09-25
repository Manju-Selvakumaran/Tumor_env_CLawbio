# Tumor_env_CLawbio

Classifying the **tumor immune microenvironment** as **hot / cold / immune-excluded**
from TCGA expression data — a ClawBio hackathon project.

Each tumor is scored on two immunology layers pulled directly from the
**cBioPortal public REST API** (targeted per-gene / per-sample queries — no bulk
dataset download), then combined into a single immune-phenotype call.

## The two layers

| Layer | Question | Genes |
|-------|----------|-------|
| **1 — Immune infiltration** | Are immune cells present? | CD8A, CD4, CD3E, FOXP3, CD68, CD163, NCAM1, CD19, CD274, PDCD1 |
| **2 — Immunosuppression** | Is the response being shut down? | TIGIT, LAG3, HAVCR2, PDCD1, CD274, IL10, TGFB1, CTLA4 |

Each layer score is the mean z-score of its panel. Decision rule (z thresholds ±0.5):

- **L1 ≥ +0.5 → HOT** (inflamed; "inflamed but immunosuppressed" if L2 ≥ +0.5)
- **L1 ≤ −0.5 → COLD** (immune desert)
- otherwise → **IMMUNE-EXCLUDED / altered**

## Two scoring scales

- **`cohort`** — within-cancer z (each tumor vs. its own cancer type). Best for a single sample.
- **`common`** — absolute batch-normalized RSEM → log2 → pooled z across **all provided
  samples**, so different cancers are directly comparable. Best for cross-cancer comparison.

## Cohorts analysed

15 samples each from **SKCM** (melanoma), **CESC** (cervical), **PAAD** (pancreatic)
TCGA PanCancer Atlas studies — 45 tumors total, validated for existence, sample type
(primary vs metastatic) and gene coverage.

**Headline result (common cross-cancer scale):** SKCM is *bimodal* (spans the coldest and
some of the hottest tumors), CESC has the largest hot fraction, and PAAD is uniformly
"lukewarm / immune-excluded" — rarely cold, rarely strongly hot.

## Repository layout

```
skills/tumor-hot-cold/tumor-hot-cold/   # packaged ClawBio skill (v0.2.0)
  ├─ SKILL.md                           #   skill definition (13/13 conformance)
  ├─ tumor_hot_cold.py                  #   live cBioPortal fetch + classifier
  ├─ data/demo_expr.json                #   offline demo cache (15 primary tumors)
  ├─ tests/                             #   pytest suite (4/4 passing)
  └─ output_live/                       #   example live run output
figures/
  ├─ hotcold_scatter.png                # L1 x L2 immune-phenotype map (45 tumors)
  └─ hotcold_heatmap.png                # gene x tumor heatmap, grouped by cancer
hotcold12.py                            # 45-sample cross-cancer analysis
hotcold12_result.json                   # full results (per-gene z, layer scores, spread)
hotcold_validation.json                 # sample existence / type / coverage check
make_figures.py                         # regenerates the figures
*_classified.csv, *_zscores.csv         # per-cohort classification tables
```

## Quick start

```bash
pip install clawbio

# Run the packaged skill (offline demo, no network)
python skills/tumor-hot-cold/tumor-hot-cold/tumor_hot_cold.py --demo --output demo/

# Classify your own samples (cross-cancer, primary tumors only)
python skills/tumor-hot-cold/tumor-hot-cold/tumor_hot_cold.py \
  --input samples.txt --output out/ --scale common --primary-only

# Reproduce the 45-tumor analysis and figures
python hotcold12.py && python make_figures.py
```

`samples.txt` is one tumor per line: `study_id<TAB or ,>sample_id`.

## References

- **Thorsson V. et al. (2018).** The Immune Landscape of Cancer. *Immunity* 48(4):812–830.e14.
  doi:10.1016/j.immuni.2018.03.023 — primary framework reference.
- Aran D., Hu Z., Butte A.J. (2017). xCell. *Genome Biology* 18:220.
- cBioPortal — https://www.cbioportal.org/
- Ayers M. et al. (2017), *J Clin Invest*; Galon J. & Bruni D. (2019), *Nat Rev Drug Discov*.

## Disclaimer

Research and educational use only. Not a medical device. Uses a transparent marker-panel
proxy for the immune-landscape framework — it does not (yet) reproduce Thorsson immune
subtypes (C1–C6) or xCell deconvolution scores (see `SKILL.md` roadmap).
