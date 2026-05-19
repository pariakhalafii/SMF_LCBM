#!/usr/bin/env bash
# 01b_merge.sh -- Merge overlapping R1/R2 with FLASH after trimming.
#
# Used when paired-end reads overlap in the middle (e.g. 2x300bp Illumina on a
# ~550bp bisulfite SMF insert -> ~50bp overlap). bwa pemerge is unreliable on
# bisulfite-converted reads; FLASH is the workhorse here.
#
# Usage: bash pipeline/01b_merge.sh config/config.example.yaml
#
# Requires: flash, python3.

set -euo pipefail

CONFIG="${1:-config/config.example.yaml}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(python3 "$HERE/_load_config.py" "$CONFIG")"

: "${sample_id:?missing sample_id}"
: "${out_dir:?missing out_dir}"
: "${threads:=4}"
: "${flash_min_overlap:=30}"
: "${flash_max_overlap:=100}"
: "${flash_max_mismatch_density:=0.25}"

TRIM_DIR="$out_dir/trimmed"
MERGE_DIR="$out_dir/merged"
mkdir -p "$MERGE_DIR"

R1_TRIM=$(ls "$TRIM_DIR"/*_val_1.fq.gz | head -n1)
R2_TRIM=$(ls "$TRIM_DIR"/*_val_2.fq.gz | head -n1)

echo "[01b] FLASH merging $R1_TRIM + $R2_TRIM"
echo "       min_overlap=$flash_min_overlap max_overlap=$flash_max_overlap mismatch=$flash_max_mismatch_density"

flash \
    --min-overlap="$flash_min_overlap" \
    --max-overlap="$flash_max_overlap" \
    --max-mismatch-density="$flash_max_mismatch_density" \
    --threads="$threads" \
    --compress \
    --output-prefix="${sample_id}.merged" \
    --output-directory="$MERGE_DIR" \
    "$R1_TRIM" "$R2_TRIM"

# FLASH writes:
#   ${sample_id}.merged.extendedFrags.fastq.gz   <- successfully merged single-end reads
#   ${sample_id}.merged.notCombined_1.fastq.gz   <- pairs that didn't merge (kept just in case)
#   ${sample_id}.merged.notCombined_2.fastq.gz
#   ${sample_id}.merged.hist / .histogram        <- distribution of merged-fragment lengths

MERGED="$MERGE_DIR/${sample_id}.merged.extendedFrags.fastq.gz"
HIST="$MERGE_DIR/${sample_id}.merged.histogram"

if [[ -f "$HIST" ]]; then
    echo "[01b] Merge length distribution (head):"
    head -n 20 "$HIST"
fi

N_MERGED=$(zcat "$MERGED" | awk 'NR%4==1' | wc -l || true)
echo "[01b] Done. Merged reads: $N_MERGED"
echo "[01b] Output: $MERGED  (single-end input for 02_align.sh)"
