# Hot vs Cold Tumor Classification (Layers 1-2)

Source: cBioPortal public REST API, `*_rna_seq_v2_mrna_median_all_sample_Zscores` (within-cohort z).
RPPA CD8/PD-L1/PD-1/pSTAT3 not in the TCGA PanCan RPPA panel; RNA covers both layers.

| Cohort | Sample | Layer 1 (infiltration) | Layer 2 (suppression) | Classification |
|---|---|---|---|---|
| SKCM | TCGA-BF-A1PU-01 | -1.131 (LOW) | -1.787 (LOW) | **COLD - immune desert** |
| CESC | TCGA-2W-A8YY-01 | -0.168 (MID) | -0.455 (MID) | **IMMUNE-EXCLUDED / altered** |
| PAAD | TCGA-2J-AAB4-01 | 0.112 (MID) | 0.081 (MID) | **IMMUNE-EXCLUDED / altered** |

Thresholds: HIGH >= 0.5, LOW <= -0.5 (z-score, within cohort).

Layer 1 panel: CD8A, CD4, FOXP3, CD68, NCAM1, CD19

Layer 2 panel: TIGIT, LAG3, HAVCR2, PDCD1, CD274, IL10, TGFB1, CTLA4