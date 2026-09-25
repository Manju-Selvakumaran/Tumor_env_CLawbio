#!/usr/bin/env python
"""
Hot vs Cold tumor classifier (Layers 1-2) from cBioPortal RNA z-scores.

Reads hotcold_raw.json (extracted from the cBioPortal public REST API:
  molecular profile *_rna_seq_v2_mrna_median_all_sample_Zscores)
and classifies each tumor by:
  Layer 1 - Immune Infiltration   (is immune infiltrate present?)
  Layer 2 - Immunosuppression      (is the response actively suppressed?)
then combines them into a hot / cold / excluded call.

z-scores are WITHIN-COHORT (per cancer type), so each score means
"relative to other tumors of the same cancer type".
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
RAW  = os.path.join(HERE, "hotcold_raw.json")

# --- panels -------------------------------------------------------------
LAYER1 = ["CD8A", "CD4", "FOXP3", "CD68", "NCAM1", "CD19"]          # infiltration
LAYER2 = ["TIGIT", "LAG3", "HAVCR2", "PDCD1", "CD274",             # suppression /
          "IL10", "TGFB1", "CTLA4"]                                # checkpoint-exhaustion
HI, LO = 0.5, -0.5                                                 # z thresholds

def band(x, hi=HI, lo=LO):
    return "HIGH" if x >= hi else ("LOW" if x <= lo else "MID")

def mean(vals):
    v = [x for x in vals if x is not None]
    return sum(v) / len(v) if v else None

def classify(l1, l2):
    """Combine the two layer scores into a phenotype."""
    if l1 >= HI:
        return ("HOT - inflamed" if l2 < HI
                else "HOT - inflamed but immunosuppressed")
    if l1 <= LO:
        return "COLD - immune desert"
    # intermediate infiltration
    return ("IMMUNE-EXCLUDED / altered (suppressed)" if l2 >= HI
            else "IMMUNE-EXCLUDED / altered")

def main():
    data = json.load(open(RAW))
    rows, report = [], []
    for tag, rec in data.items():
        rna = rec["RNA"]
        l1v = [rna.get(g) for g in LAYER1]
        l2v = [rna.get(g) for g in LAYER2]
        l1, l2 = mean(l1v), mean(l2v)
        call = classify(l1, l2)
        rows.append((tag, rec["sample"], l1, l2, call))
        report.append({
            "cohort": tag, "sample": rec["sample"],
            "layer1_infiltration_score": round(l1, 3), "layer1_band": band(l1),
            "layer2_suppression_score":  round(l2, 3), "layer2_band": band(l2),
            "layer1_markers": {g: rna.get(g) for g in LAYER1},
            "layer2_markers": {g: rna.get(g) for g in LAYER2},
            "classification": call,
        })

    # console table
    print(f"{'Cohort':6} {'Sample':17} {'L1 infil':>9} {'L2 suppr':>9}  Classification")
    print("-" * 78)
    for tag, samp, l1, l2, call in rows:
        print(f"{tag:6} {samp:17} {l1:>9.2f} {l2:>9.2f}  {call}")

    json.dump(report, open(os.path.join(HERE, "hotcold_result.json"), "w"), indent=2)

    # markdown
    md = ["# Hot vs Cold Tumor Classification (Layers 1-2)",
          "",
          "Source: cBioPortal public REST API, "
          "`*_rna_seq_v2_mrna_median_all_sample_Zscores` (within-cohort z).",
          "RPPA CD8/PD-L1/PD-1/pSTAT3 not in the TCGA PanCan RPPA panel; "
          "RNA covers both layers.",
          "",
          "| Cohort | Sample | Layer 1 (infiltration) | Layer 2 (suppression) | Classification |",
          "|---|---|---|---|---|"]
    for r in report:
        md.append(f"| {r['cohort']} | {r['sample']} | "
                  f"{r['layer1_infiltration_score']} ({r['layer1_band']}) | "
                  f"{r['layer2_suppression_score']} ({r['layer2_band']}) | "
                  f"**{r['classification']}** |")
    md += ["",
           f"Thresholds: HIGH >= {HI}, LOW <= {LO} (z-score, within cohort).",
           "",
           "Layer 1 panel: " + ", ".join(LAYER1),
           "",
           "Layer 2 panel: " + ", ".join(LAYER2)]
    open(os.path.join(HERE, "hotcold_report.md"), "w").write("\n".join(md))
    print("\nWrote hotcold_result.json and hotcold_report.md")

if __name__ == "__main__":
    main()
