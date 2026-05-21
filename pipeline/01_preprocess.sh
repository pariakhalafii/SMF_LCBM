#!/usr/bin/env bash
# 01_preprocess.sh -- FastQC + Trim Galore for paired-end SMF reads.
#
# Usage: bash pipeline/01_preprocess.sh config/config.example.yaml
#
# Requires: fastqc, trim_galore (which itself requires cutadapt), python3.

set -euo pipefail

CONFIG="${1:-config/config.example.yaml}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(python3 "$HERE/_load_config.py" "$CONFIG")"

: "${sample_id:?missing sample_id}"
: "${fastq_r1:?missing fastq_r1}"
: "${fastq_r2:?missing fastq_r2}"
: "${out_dir:?missing out_dir}"
: "${threads:=4}"

QC_DIR="$out_dir/qc"
TRIM_DIR="$out_dir/trimmed"
mkdir -p "$QC_DIR" "$TRIM_DIR"

echo "[01] FastQC on raw reads"
fastqc -t "$threads" -o "$QC_DIR" "$fastq_r1" "$fastq_r2"

echo "[01] Trim Galore (paired-end, bisulfite-aware adapter trimming)"
trim_galore \
    --paired \
    --cores "$threads" \
    --output_dir "$TRIM_DIR" \
    --fastqc \
    --length 70 \
    --stringency 6 \
    "$fastq_r1" "$fastq_r2"

echo "[01] Done. Trimmed FASTQs in $TRIM_DIR"
