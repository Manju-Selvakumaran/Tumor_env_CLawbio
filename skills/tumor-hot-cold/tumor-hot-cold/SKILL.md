---
name: tumor-hot-cold
description: >-
  Classify TCGA tumors as hot / cold / immune-excluded from cBioPortal RNA expression z-scores across two immunology layers (immune infiltration and immunosuppression), with no bulk dataset download.
license: MIT
metadata:
  version: "0.2.0"
  author: "ClawBio Hackathon"
  domain: cancer-immunology
  tags:
    - cancer
    - immuno-oncology
    - TCGA
    - cBioPortal
    - immune-phenotype
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🦖"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: python
    trigger_keywords:
        - hot cold tumor
        - hot vs cold
        - immune infiltration classification
        - immune excluded
        - immunosuppression score
        - tumor immune phenotype
        - cBioPortal z-score classification
---

# 🦖 Tumor Hot Cold

You are **Tumor Hot Cold**, a specialised ClawBio agent for cancer-immunology. Your role is to Classify TCGA tumors as hot / cold / immune-excluded from cBioPortal RNA expression z-scores across two immunology layers (immune infiltration and immunosuppression), with no bulk dataset download.

## Why This Exists

- **Without it**: A researcher must open cBioPortal, query each gene panel one sample at a time, copy z-score tables out by hand, then eyeball whether a tumor is "hot" or "cold" with no consistent rule — slow, error-prone, and irreproducible.
- **With it**: Give a study + sample id and get a per-layer immune score and a hot / cold / immune-excluded call in seconds, with a report and reproducibility bundle — no bulk dataset download.
- **Why ClawBio**: Scores come from real cBioPortal expression z-scores via the public REST API and a fixed, transparent marker panel and threshold — not an LLM guessing from gene names.

## Core Capabilities

1. **Targeted data fetch**: Pull per-gene, per-sample mRNA expression from the cBioPortal public REST API (targeted queries, no bulk download).
2. **Layer 1 — immune infiltration** (10 genes): Score the infiltration panel (CD8A, CD4, CD3E, FOXP3, CD68, CD163, NCAM1, CD19, CD274, PDCD1) — is immune infiltrate present?
3. **Layer 2 — immunosuppression**: Score the suppression / checkpoint-exhaustion panel (TIGIT, LAG3, HAVCR2, PDCD1, CD274, IL10, TGFB1, CTLA4) — is the response being actively suppressed?
4. **Two scoring scales**: `cohort` (within-cancer z, best for a single sample) or `common` (absolute batch-normalized RSEM → log2 → pooled z across all provided samples, so different cancers are directly comparable). `auto` picks `common` when ≥3 samples.
5. **Validation preflight**: Verify every sample exists, report its type (primary / metastatic / normal), and check per-gene data coverage; `--primary-only` drops non-primary samples for apples-to-apples comparison.
6. **Cross-cancer comparison**: Per-cancer spread summary (mean / min / max / SD and hot·excluded·cold counts) when multiple samples are supplied.
7. **Combined phenotype + reproducible outputs**: Fold both layers into a hot / cold / immune-excluded call per tumor; emit `report.md`, `result.json` and a reproducibility bundle.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| cBioPortal study + sample identifiers | `.txt` | study_id, sample_id (one per line, tab or comma separated) | `demo_samples.txt` |

## Workflow

When the user asks to classify a tumor's immune phenotype (hot vs cold) from TCGA / cBioPortal expression data:

1. **Validate**: Read the samples file (`study_id` + `sample_id` per line); preflight each sample for existence, type, and gene coverage; optionally drop non-primary samples.
2. **Fetch**: `common` scale → absolute RSEM (`*_rna_seq_v2_mrna`), log2, pooled-z across all provided samples; `cohort` scale → within-cancer z (`*_rna_seq_v2_mrna_median_all_sample_Zscores`).
3. **Score & classify**: Average each panel to a Layer 1 / Layer 2 score and apply the decision rule; build a per-cancer spread summary.
4. **Report**: Write `result.json` (validation, per-marker z, layer scores, bands, classification, per-cancer summary) and `report.md`, plus the reproducibility bundle (`commands.sh`, `environment.yml`, `checksums.sha256`).

## CLI Reference

```bash
# Standard usage (auto scale)
python skills/tumor-hot-cold/tumor_hot_cold.py \
  --input <input_file> --output <report_dir>

# Cross-cancer comparison on a common scale, primary tumors only
python skills/tumor-hot-cold/tumor_hot_cold.py \
  --input <input_file> --output <report_dir> --scale common --primary-only

# Single sample relative to its own cancer type
python skills/tumor-hot-cold/tumor_hot_cold.py \
  --input <input_file> --output <report_dir> --scale cohort

# Demo mode (bundled cache, no network)
python skills/tumor-hot-cold/tumor_hot_cold.py --demo --output /tmp/tumor-hot-cold_demo

# Via ClawBio runner
python clawbio.py run hotcold --input <file> --output <dir>
python clawbio.py run hotcold --demo
```

## Demo

```bash
python clawbio.py run hotcold --demo
```

Expected output: a fully offline run (bundled expression cache, no network) classifying **15 primary TCGA tumors** (5 each from SKCM, CESC, PAAD) on the common cross-cancer scale — e.g. SKCM `TCGA-BF-A1PU-01` → **COLD (immune desert)**, SKCM `TCGA-BF-A1PX-01` → **HOT (inflamed but immunosuppressed)**, PAAD `TCGA-2J-AABP-01` → **HOT** — plus a per-cancer spread table, `report.md`, `result.json` and a reproducibility bundle. Demonstrates SKCM's bimodal hot/cold spread and PAAD's tighter distribution.

## Algorithm / Methodology

Describe the core methodology so an AI agent can apply it even without the Python script:

1. **Preflight**: Map panel HGNC symbols → Entrez (`/genes/fetch`); verify each sample via `/samples/fetch` (existence + `sampleType`); optionally keep only `Primary Solid Tumor`.
2. **Fetch (scale-dependent)**: `common` → absolute RSEM from `*_rna_seq_v2_mrna`, transform `log2(x+1)`, then standardize each gene across **all provided samples** (pooled z) so cancers share one scale. `cohort` → pull within-cancer z directly from `*_rna_seq_v2_mrna_median_all_sample_Zscores`.
3. **Score layers**: Layer 1 (infiltration) = mean z of {CD8A, CD4, CD3E, FOXP3, CD68, CD163, NCAM1, CD19, CD274, PDCD1}; Layer 2 (suppression) = mean z of {TIGIT, LAG3, HAVCR2, PDCD1, CD274, IL10, TGFB1, CTLA4}.
4. **Classify**: L1 ≥ +0.5 → HOT (inflamed; "inflamed but immunosuppressed" if L2 ≥ +0.5); L1 ≤ −0.5 → COLD (immune desert); otherwise IMMUNE-EXCLUDED / altered. Summarize spread per cancer.

**Key thresholds / parameters**:
- HIGH band: z ≥ +0.5; LOW band: z ≤ −0.5.
- `cohort` z is relative to a tumor's own cancer type; `common` pooled-z is relative to the set of samples supplied in the run (enables cross-cancer comparison).
- Marker panels and the hot/cold immune-phenotype framework follow the TCGA immune-landscape literature (Thorsson et al. 2018; Galon immune contexture; Ayers IFN-γ inflamed signature).

**Roadmap / limitations**: This release scores a transparent marker-panel proxy. It does **not** yet reproduce Thorsson's six immune subtypes (C1–C6) or xCell deconvolution scores — those are not exposed by the cBioPortal REST API for these studies (they require the Thorsson supplementary table / running xCell). A planned extension cross-references each call against Thorsson subtype + xCell ImmuneScore so each label rests on two independent sources.

## Example Queries

- "Is TCGA-BF-A1PU-01 in the SKCM cohort a hot or cold tumor?"
- "Classify the immune phenotype of these cBioPortal samples across infiltration and immunosuppression layers."

## Output Structure

```
output_directory/
├── report.md              # Per-tumor table + per-cancer spread + validation summary
├── result.json            # Validation, per-marker z, layer scores, classification, per-cancer summary
└── reproducibility/
    ├── commands.sh        # Exact commands to reproduce
    ├── environment.yml    # Conda/pip environment snapshot
    └── checksums.sha256   # SHA-256 of outputs, relative to output_directory
```

Verify the run with
`cd <output_directory> && sha256sum -c reproducibility/checksums.sha256`.

## Dependencies

**Required** (in `requirements.txt` or skill-level install):
- `python >= 3.11` — stdlib only (`urllib`, `json`, `argparse`); the cBioPortal REST calls and scoring need no third-party packages.

## Safety

- **Local-first**: No data upload without explicit consent
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: Log all operations to reproducibility bundle
- **No hallucinated science**: All parameters trace to cited databases

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- hot cold tumor
- hot vs cold
- immune infiltration classification

**Chaining partners** — this skill connects with:
- `xena-tcga-gene-query`: supplies additional TCGA expression/survival context (e.g. methylation or extra markers) for samples classified here.
- `omics-target-evidence-mapper`: takes the hot/cold call as immune-context input when triaging immuno-oncology targets.

## Citations

- **Thorsson V. et al. (2018), *Immunity* 48(4):812–830.e14** — "The Immune Landscape of Cancer." doi:10.1016/j.immuni.2018.03.023 — primary reference for the TCGA tumor immune-landscape framework this skill operationalizes.
- Aran D., Hu Z. & Butte A.J. (2017), *Genome Biology* 18:220 — **xCell**: cell-type enrichment from bulk transcriptomes (planned second-source ImmuneScore).
- [cBioPortal](https://www.cbioportal.org/) — public REST API and TCGA PanCancer Atlas mRNA expression profiles (absolute RSEM and z-score).
- Ayers M. et al. (2017), *J Clin Invest* — IFN-γ–related mRNA profile defining the T-cell-inflamed ("hot") tumor phenotype.
- Galon J. & Bruni D. (2019), *Nat Rev Drug Discov* — immune contexture and hot/cold/excluded tumor classification framework.
