# Tumor Hot vs Cold Classification (Layers 1-2)

**Generated**: 2026-09-25 13:26:27
**Scoring scale**: `common`
**Source**: cBioPortal live: pooled z across provided samples from log2 batch-normalized RSEM (_rna_seq_v2_mrna)

## Validation preflight
- Samples requested: 9 | scored: 9
- Sample types: {'Primary Solid Tumor': 9}
- Gene coverage: 9/9 samples have all 16 genes
- Dropped: none

## Per-tumor classification

| Cancer | Sample | Type | L1 (infiltration) | L2 (suppression) | Classification |
|---|---|---|---|---|---|
| SKCM | TCGA-BF-A1PU-01 | Primary Solid Tumor | -1.326 (LOW) | -1.719 (LOW) | **COLD - immune desert** |
| SKCM | TCGA-BF-A1PX-01 | Primary Solid Tumor | 0.748 (HIGH) | 0.841 (HIGH) | **HOT - inflamed but immunosuppressed** |
| SKCM | TCGA-BF-A1Q0-01 | Primary Solid Tumor | 0.359 (MID) | -0.029 (MID) | **IMMUNE-EXCLUDED / altered** |
| CESC | TCGA-C5-A1BK-01 | Primary Solid Tumor | 0.678 (HIGH) | 0.811 (HIGH) | **HOT - inflamed but immunosuppressed** |
| CESC | TCGA-C5-A1BN-01 | Primary Solid Tumor | -1.473 (LOW) | -0.993 (LOW) | **COLD - immune desert** |
| CESC | TCGA-2W-A8YY-01 | Primary Solid Tumor | -0.342 (MID) | -0.32 (MID) | **IMMUNE-EXCLUDED / altered** |
| PAAD | TCGA-2J-AAB4-01 | Primary Solid Tumor | 0.056 (MID) | -0.066 (MID) | **IMMUNE-EXCLUDED / altered** |
| PAAD | TCGA-2J-AABP-01 | Primary Solid Tumor | 0.713 (HIGH) | 0.999 (HIGH) | **HOT - inflamed but immunosuppressed** |
| PAAD | TCGA-2J-AABR-01 | Primary Solid Tumor | 0.587 (HIGH) | 0.475 (MID) | **HOT - inflamed** |

## Per-cancer spread & cross-cancer comparison

| Cancer | n | L1 mean | L1 min | L1 max | L1 SD | hot | excl | cold |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| PAAD | 3 | 0.452 | 0.056 | 0.713 | 0.285 | 2 | 1 | 0 |
| SKCM | 3 | -0.073 | -1.326 | 0.748 | 0.9 | 1 | 1 | 1 |
| CESC | 3 | -0.379 | -1.473 | 0.678 | 0.879 | 1 | 1 | 1 |

Thresholds: HIGH >= 0.5, LOW <= -0.5.

**Layer 1 — immune infiltration**: CD8A, CD4, CD3E, FOXP3, CD68, CD163, NCAM1, CD19, CD274, PDCD1

**Layer 2 — immunosuppression**: TIGIT, LAG3, HAVCR2, PDCD1, CD274, IL10, TGFB1, CTLA4

Marker panels and the hot/cold immune-phenotype framework follow the TCGA immune-landscape literature (Thorsson et al., *Immunity* 2018); this skill uses a transparent marker-panel proxy rather than the full immune-subtype model or xCell deconvolution (see SKILL.md roadmap).

---

> **Disclaimer**: Research and educational use only. Not a medical device.