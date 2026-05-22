#!/usr/bin/env bash
# se_01_preprocess.sh -- FastQC + Trim Galore for SINGLE-END (merged) SMF reads.
#
# Usage: bash pipeline/se_01_preprocess.sh config/config.example1.yaml
#
# Requires: fastqc, trim_galore (which itself requires cutadapt), python3.

set -euo pipefail

CONFIG="${1:-config/config.example1.yaml}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(python3 "$HERE/_load_config.py" "$CONFIG")"

: "${sample_id:?missing sample_id}"
: "${fastq_se:?missing fastq_se}"
: "${out_dir:?missing out_dir}"
: "${threads:=4}"

QC_DIR="$out_dir/qc"
TRIM_DIR="$out_dir/trimmed"
mkdir -p "$QC_DIR" "$TRIM_DIR"

echo "[se_01] FastQC on raw reads"
fastqc -t "$threads" -o "$QC_DIR" "$fastq_se"

echo "[se_01] Trim Galore (single-end, bisulfite-aware adapter trimming)"
trim_galore \
    --cores "$threads" \
    --output_dir "$TRIM_DIR" \
    --fastqc \
    --stringency 6 \
    "$fastq_se"

echo "[se_01] Done. Trimmed FASTQ in $TRIM_DIR"
