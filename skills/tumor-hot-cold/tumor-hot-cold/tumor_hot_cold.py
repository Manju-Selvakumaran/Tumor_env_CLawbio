#!/usr/bin/env python3
"""
tumor_hot_cold.py — ClawBio Tumor Hot Cold Skill
==================================================
Classify TCGA tumors as hot / cold / immune-excluded from cBioPortal RNA
expression across two immunology layers (immune infiltration and
immunosuppression), with no bulk dataset download.

Two scoring scales:
  * cohort  — each sample's within-cohort z (relative to its own cancer type;
              profile *_rna_seq_v2_mrna_median_all_sample_Zscores). Best for a
              single sample.
  * common  — absolute batch-normalized RSEM -> log2 -> pooled z across ALL the
              provided samples, so different cancers are directly comparable.
              Best for multi-sample / cross-cancer comparison.
  (default: auto -> common when >=3 samples, else cohort)

Every run performs a validation preflight: sample existence, sample type
(primary / metastatic / normal) and per-gene data coverage. Use --primary-only
to keep just Primary Solid Tumor samples for an apples-to-apples comparison.

Author:  ClawBio Hackathon
Version: 0.2.0

Usage:
    python tumor_hot_cold.py --input samples.txt --output out/
    python tumor_hot_cold.py --input samples.txt --output out/ --scale common --primary-only
    python tumor_hot_cold.py --demo --output demo/

Input file: one tumor per line, "study_id<TAB or ,>sample_id"; '#' comments ignored.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics as st
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

# ---------------------------------------------------------------------------
# Panels & configuration
# ---------------------------------------------------------------------------
# Layer 1 — immune infiltration (T / Treg / macrophage / NK / B + activity)
LAYER1 = ["CD8A", "CD4", "CD3E", "FOXP3", "CD68", "CD163",
          "NCAM1", "CD19", "CD274", "PDCD1"]
# Layer 2 — immunosuppression / checkpoint-exhaustion
LAYER2 = ["TIGIT", "LAG3", "HAVCR2", "PDCD1", "CD274", "IL10", "TGFB1", "CTLA4"]
GENES = sorted(set(LAYER1) | set(LAYER2))
HI, LO = 0.5, -0.5

CBIO = "https://www.cbioportal.org/api"
ABS_PROFILE = "_rna_seq_v2_mrna"                              # batch-normalized RSEM
COHORTZ_PROFILE = "_rna_seq_v2_mrna_median_all_sample_Zscores"  # within-cohort z
DEMO_CACHE = SKILL_DIR / "data" / "demo_expr.json"

# ---------------------------------------------------------------------------
# cBioPortal REST helpers (targeted queries — no bulk download)
# ---------------------------------------------------------------------------

def _post(path: str, body):
    req = urllib.request.Request(
        CBIO + path, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=60))


def _gene_maps():
    g = _post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", GENES)
    sym2ez = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in g}
    return sym2ez, {v: k for k, v in sym2ez.items()}


def fetch_sample_meta(study: str, samps: list) -> dict:
    """Return {sample: {'sample_type':..., 'exists':bool}} via samples/fetch."""
    found = _post("/samples/fetch?projection=DETAILED",
                  {"sampleIdentifiers": [{"sampleId": s, "studyId": study} for s in samps]})
    fmap = {x["sampleId"]: x for x in found}
    return {s: {"sample_type": fmap.get(s, {}).get("sampleType"),
                "exists": s in fmap} for s in samps}


def fetch_abs_log2(study: str, samps: list, sym2ez, ez2sym) -> dict:
    """{sample: {gene: log2(RSEM+1)}} from the absolute profile."""
    recs = _post(f"/molecular-profiles/{study}{ABS_PROFILE}/molecular-data/fetch?projection=SUMMARY",
                 {"sampleIds": samps, "entrezGeneIds": list(sym2ez.values())})
    out = {s: {} for s in samps}
    for r in recs:
        gene = ez2sym.get(r["entrezGeneId"]); v = r.get("value")
        if gene and v is not None and r["sampleId"] in out:
            out[r["sampleId"]][gene] = math.log2(v + 1.0)
    return out


def fetch_cohort_z(study: str, samps: list, sym2ez, ez2sym) -> dict:
    """{sample: {gene: within-cohort z}} from the all-sample z-score profile."""
    recs = _post(f"/molecular-profiles/{study}{COHORTZ_PROFILE}/molecular-data/fetch?projection=SUMMARY",
                 {"sampleIds": samps, "entrezGeneIds": list(sym2ez.values())})
    out = {s: {} for s in samps}
    for r in recs:
        gene = ez2sym.get(r["entrezGeneId"]); v = r.get("value")
        if gene and v is not None and r["sampleId"] in out:
            out[r["sampleId"]][gene] = v
    return out

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def pooled_z(expr_by_key: dict) -> dict:
    """Standardize each gene across ALL provided samples (common scale)."""
    keys = list(expr_by_key)
    stats = {}
    for gene in GENES:
        vals = [expr_by_key[k][gene] for k in keys if gene in expr_by_key[k]]
        if vals:
            stats[gene] = (st.mean(vals), st.pstdev(vals) or 1e-9)
    return {k: {g: (expr_by_key[k][g] - stats[g][0]) / stats[g][1]
                for g in expr_by_key[k] if g in stats} for k in keys}


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def _band(x):
    return "HIGH" if x >= HI else ("LOW" if x <= LO else "MID")


def classify(l1, l2):
    if l1 >= HI:
        return "HOT - inflamed" if l2 < HI else "HOT - inflamed but immunosuppressed"
    if l1 <= LO:
        return "COLD - immune desert"
    return ("IMMUNE-EXCLUDED / altered (suppressed)" if l2 >= HI
            else "IMMUNE-EXCLUDED / altered")


def build_records(z_by_key: dict, meta: dict) -> list:
    recs = []
    for key, zr in z_by_key.items():
        cancer, sample = key.split("|", 1)
        l1 = _mean(zr.get(g) for g in LAYER1)
        l2 = _mean(zr.get(g) for g in LAYER2)
        recs.append({
            "cancer": cancer, "sample": sample,
            "sample_type": meta.get(key, {}).get("sample_type"),
            "layer1_infiltration_score": round(l1, 3), "layer1_band": _band(l1),
            "layer2_suppression_score": round(l2, 3), "layer2_band": _band(l2),
            "genes_used": len(zr),
            "z": {g: round(zr[g], 3) for g in GENES if g in zr},
            "classification": classify(l1, l2),
        })
    return recs


def per_cancer_summary(recs: list) -> dict:
    out = {}
    for lab in sorted({r["cancer"] for r in recs}):
        l1 = [r["layer1_infiltration_score"] for r in recs if r["cancer"] == lab]
        l2 = [r["layer2_suppression_score"] for r in recs if r["cancer"] == lab]
        calls = [r["classification"] for r in recs if r["cancer"] == lab]
        out[lab] = {
            "n": len(l1),
            "L1_mean": round(st.mean(l1), 3), "L1_min": min(l1), "L1_max": max(l1),
            "L1_sd": round(st.pstdev(l1), 3),
            "L2_mean": round(st.mean(l2), 3), "L2_min": min(l2), "L2_max": max(l2),
            "L2_sd": round(st.pstdev(l2), 3),
            "hot": sum("HOT" in c for c in calls),
            "cold": sum("COLD" in c for c in calls),
            "excluded": sum("EXCLUDED" in c for c in calls),
        }
    return out

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def parse_input(input_path: Path):
    pairs = []
    for raw in input_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.replace(",", "\t").split("\t") if p.strip()]
        if len(parts) >= 2:
            pairs.append((parts[0], parts[1]))
    return pairs


def _write_outputs(output_dir, records, summary, validation, scale,
                   invocation, source_note) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    multi = len(records) > 1
    results = {
        "skill": "tumor-hot-cold",
        "generated_at": datetime.now().isoformat(),
        "scale": scale, "source": source_note,
        "layer1_panel": LAYER1, "layer2_panel": LAYER2,
        "thresholds": {"high": HI, "low": LO},
        "validation": validation,
        "tumors": records,
        "per_cancer": summary if multi else {},
    }
    (output_dir / "result.json").write_text(
        json.dumps(results, indent=2, default=str), encoding="utf-8")

    md = [
        "# Tumor Hot vs Cold Classification (Layers 1-2)",
        "",
        f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Scoring scale**: `{scale}`",
        f"**Source**: {source_note}",
        "",
        "## Validation preflight",
        f"- Samples requested: {validation['n_requested']} | "
        f"scored: {validation['n_scored']}",
        f"- Sample types: {validation['sample_types']}",
        f"- Gene coverage: {validation['gene_coverage']}",
        f"- Dropped: {validation['dropped'] or 'none'}",
        "",
        "## Per-tumor classification",
        "",
        "| Cancer | Sample | Type | L1 (infiltration) | L2 (suppression) | Classification |",
        "|---|---|---|---|---|---|",
    ]
    for r in records:
        md.append(
            f"| {r['cancer']} | {r['sample']} | {r['sample_type'] or '?'} | "
            f"{r['layer1_infiltration_score']} ({r['layer1_band']}) | "
            f"{r['layer2_suppression_score']} ({r['layer2_band']}) | "
            f"**{r['classification']}** |")
    if len(records) > 1:
        md += ["", "## Per-cancer spread & cross-cancer comparison", "",
               "| Cancer | n | L1 mean | L1 min | L1 max | L1 SD | hot | excl | cold |",
               "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
        for lab, s in sorted(summary.items(), key=lambda kv: -kv[1]["L1_mean"]):
            md.append(f"| {lab} | {s['n']} | {s['L1_mean']} | {s['L1_min']} | "
                      f"{s['L1_max']} | {s['L1_sd']} | {s['hot']} | "
                      f"{s['excluded']} | {s['cold']} |")
    md += [
        "",
        f"Thresholds: HIGH >= {HI}, LOW <= {LO}.",
        "",
        f"**Layer 1 — immune infiltration**: {', '.join(LAYER1)}",
        "",
        f"**Layer 2 — immunosuppression**: {', '.join(LAYER2)}",
        "",
        "Marker panels and the hot/cold immune-phenotype framework follow the "
        "TCGA immune-landscape literature (Thorsson et al., *Immunity* 2018); "
        "this skill uses a transparent marker-panel proxy rather than the full "
        "immune-subtype model or xCell deconvolution (see SKILL.md roadmap).",
        "",
        "---",
        "",
        "> **Disclaimer**: Research and educational use only. Not a medical device.",
    ]
    (output_dir / "report.md").write_text("\n".join(md), encoding="utf-8")

    write_commands_sh(
        output_dir,
        "# Tumor Hot Cold - reproducibility\n"
        "set -euo pipefail\n"
        "conda env create -f environment.yml\n"
        "conda activate tumor-hot-cold\n"
        f"python tumor_hot_cold.py {invocation}\n"
        '( cd "$(dirname "$0")/.." && sha256sum -c reproducibility/checksums.sha256 )')
    write_environment_yml(output_dir, env_name="tumor-hot-cold",
                          python_version="3.11", pip_deps=["python >= 3.11"],
                          conda_deps=[])
    write_checksums([output_dir / "report.md", output_dir / "result.json"],
                    output_dir, anchor=output_dir)
    return results

# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def run(input_path: Path, output_dir: Path, scale: str = "auto",
        primary_only: bool = False) -> dict:
    pairs = parse_input(input_path)
    if not pairs:
        raise ValueError("No 'study_id  sample_id' pairs found in input file.")
    sym2ez, ez2sym = _gene_maps()

    by_study = {}
    for study, sample in pairs:
        by_study.setdefault(study, []).append(sample)

    # --- validation preflight: existence + sample type ---
    meta, dropped = {}, []
    kept_by_study = {}
    for study, samps in by_study.items():
        lab = study.split("_")[0].upper()
        smeta = fetch_sample_meta(study, samps)
        kept = []
        for s in samps:
            key = f"{lab}|{s}"
            if not smeta[s]["exists"]:
                dropped.append(f"{key} (not found)"); continue
            if primary_only and smeta[s]["sample_type"] != "Primary Solid Tumor":
                dropped.append(f"{key} ({smeta[s]['sample_type']}, non-primary)"); continue
            meta[key] = {"study": study, "sample_type": smeta[s]["sample_type"]}
            kept.append(s)
        if kept:
            kept_by_study[study] = kept

    n_kept = sum(len(v) for v in kept_by_study.values())
    if n_kept == 0:
        raise ValueError(f"No samples left after validation. Dropped: {dropped}")

    chosen = scale
    if scale == "auto":
        chosen = "common" if n_kept >= 3 else "cohort"

    # --- fetch + score ---
    expr_by_key = {}    # for common scale (log2) / used only to check coverage in cohort
    z_by_key = {}
    for study, samps in kept_by_study.items():
        lab = study.split("_")[0].upper()
        if chosen == "common":
            log2 = fetch_abs_log2(study, samps, sym2ez, ez2sym)
            for s in samps:
                expr_by_key[f"{lab}|{s}"] = log2[s]
        else:  # cohort
            zc = fetch_cohort_z(study, samps, sym2ez, ez2sym)
            for s in samps:
                z_by_key[f"{lab}|{s}"] = zc[s]

    if chosen == "common":
        z_by_key = pooled_z(expr_by_key)
        coverage_src = expr_by_key
    else:
        coverage_src = z_by_key

    # gene coverage check
    incomplete = [k for k, d in coverage_src.items() if len(d) != len(GENES)]
    gene_coverage = (f"{sum(len(d) == len(GENES) for d in coverage_src.values())}"
                     f"/{len(coverage_src)} samples have all {len(GENES)} genes")

    records = build_records(z_by_key, meta)
    summary = per_cancer_summary(records)
    type_counts = {}
    for k in meta:
        t = meta[k]["sample_type"]; type_counts[t] = type_counts.get(t, 0) + 1
    validation = {
        "n_requested": len(pairs), "n_scored": len(records),
        "sample_types": type_counts, "gene_coverage": gene_coverage,
        "incomplete_samples": incomplete, "dropped": dropped,
        "primary_only": primary_only,
    }
    inv = (f"--input {input_path} --output {output_dir} "
           f"--scale {chosen}" + (" --primary-only" if primary_only else ""))
    src = ("cBioPortal live: pooled z across provided samples from "
           f"log2 batch-normalized RSEM ({ABS_PROFILE})" if chosen == "common"
           else f"cBioPortal live: within-cohort z ({COHORTZ_PROFILE})")
    return _write_outputs(output_dir, records, summary, validation, chosen, inv, src)

# ---------------------------------------------------------------------------
# Demo mode (offline — uses bundled data/demo_expr.json, no network)
# ---------------------------------------------------------------------------

def run_demo(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = json.loads(DEMO_CACHE.read_text(encoding="utf-8"))
    expr_by_key = cache["expr"]
    meta = {k: {"study": v["study"], "sample_type": v["sample_type"]}
            for k, v in cache["meta"].items()}
    z_by_key = pooled_z(expr_by_key)
    records = build_records(z_by_key, meta)
    summary = per_cancer_summary(records)
    type_counts = {}
    for k in meta:
        t = meta[k]["sample_type"]; type_counts[t] = type_counts.get(t, 0) + 1
    validation = {
        "n_requested": len(expr_by_key), "n_scored": len(records),
        "sample_types": type_counts,
        "gene_coverage": f"{len(records)}/{len(records)} samples have all {len(GENES)} genes",
        "incomplete_samples": [], "dropped": [], "primary_only": True,
    }
    results = _write_outputs(
        output_dir, records, summary, validation, "common",
        "--demo --output " + str(output_dir),
        "bundled offline demo cache (5 primary tumors x 3 cancers: SKCM, CESC, PAAD)")
    print(f"  Demo complete ({len(records)} tumors, common cross-cancer scale). "
          f"Output: {output_dir}")
    for r in records:
        print(f"    {r['cancer']:5} {r['sample']:17} "
              f"L1={r['layer1_infiltration_score']:>6} "
              f"L2={r['layer2_suppression_score']:>6}  -> {r['classification']}")
    return results

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify TCGA tumors hot / cold / immune-excluded from "
                    "cBioPortal RNA (Layer 1 infiltration + Layer 2 immunosuppression).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", dest="input_path",
                        help="samples file: 'study_id  sample_id' per line")
    parser.add_argument("--output", dest="output_dir",
                        default="/tmp/tumor-hot-cold_output", help="Output directory")
    parser.add_argument("--scale", choices=["auto", "cohort", "common"], default="auto",
                        help="cohort=within-cancer z; common=cross-cancer pooled z "
                             "(auto: common if >=3 samples)")
    parser.add_argument("--primary-only", action="store_true",
                        help="keep only Primary Solid Tumor samples")
    parser.add_argument("--demo", action="store_true",
                        help="Run offline with bundled demo data")
    args = parser.parse_args()
    out = Path(args.output_dir)

    if args.demo:
        run_demo(out)
    elif args.input_path:
        run(Path(args.input_path), out, scale=args.scale,
            primary_only=args.primary_only)
        print(f"  Done. Output: {out}\n  Report: {out / 'report.md'}")
    else:
        parser.error("Provide --input <file> or --demo")


if __name__ == "__main__":
    main()
