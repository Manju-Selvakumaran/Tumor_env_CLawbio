#!/usr/bin/env bash
# Tumor Hot Cold - reproducibility
set -euo pipefail
conda env create -f environment.yml
conda activate tumor-hot-cold
python tumor_hot_cold.py --input examples\demo_samples.txt --output C:\Users\manju\clawbio_hackathon\skills\tumor-hot-cold\tumor-hot-cold\output_live --scale common --primary-only
( cd "$(dirname "$0")/.." && sha256sum -c reproducibility/checksums.sha256 )
