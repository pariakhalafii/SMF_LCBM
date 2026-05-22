#!/usr/bin/env bash
# se_02_align.sh -- Bismark SINGLE-END alignment for merged SMF reads.
#
# Usage: bash pipeline/se_02_align.sh config/config.example1.yaml
#
# Requires: bismark, bowtie2, samtools, python3.

set -euo pipefail

CONFIG="${1:-config/config.example1.yaml}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(python3 "$HERE/_load_config.py" "$CONFIG")"

: "${sample_id:?missing sample_id}"
: "${genome_dir:?missing genome_dir}"
: "${out_dir:?missing out_dir}"
: "${threads:=4}"
: "${bismark_options:=--non_directional --un --ambiguous}"
: "${deduplicate:=no}"

TRIM_DIR="$out_dir/trimmed"
ALIGN_DIR="$out_dir/align"
mkdir -p "$ALIGN_DIR"

# Trim Galore single-end output keeps the original name with _trimmed suffix.
SE_TRIM=$(ls "$TRIM_DIR"/*_trimmed.fq* 2>/dev/null | head -n1)

if [[ -z "$SE_TRIM" ]]; then
    echo "ERROR: trimmed FASTQ not found in $TRIM_DIR" >&2
    echo "       Run pipeline/se_01_preprocess.sh first." >&2
    exit 1
fi

if [[ ! -f "$genome_dir/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2" ]]; then
    echo "[se_02] Bismark genome preparation (one-time)"
    bismark_genome_preparation --bowtie2 "$genome_dir"
fi

echo "[se_02] Bismark single-end alignment"
echo "       Input:   $SE_TRIM"
echo "       Options: $bismark_options"
bismark \
    --genome "$genome_dir" \
    "$SE_TRIM" \
    -p "$threads" \
    $bismark_options \
    -o "$ALIGN_DIR"

# Single-end Bismark output is *_bismark_bt2.bam
BAM=$(ls "$ALIGN_DIR"/*_bismark_bt2.bam 2>/dev/null | head -n1)

DEDUP_LC=$(printf '%s' "$deduplicate" | tr '[:upper:]' '[:lower:]')
if [[ "$DEDUP_LC" == "yes" || "$DEDUP_LC" == "true" ]]; then
    echo "[se_02] Deduplicating $BAM"
    deduplicate_bismark --single --bam --output_dir "$ALIGN_DIR" "$BAM"
    BAM=$(ls "$ALIGN_DIR"/*deduplicated.bam | head -n1)
else
    echo "[se_02] Skipping deduplication (amplicon mode)"
fi

SORTED="${BAM%.bam}.sorted.bam"
echo "[se_02] Sorting & indexing"
samtools sort -@ "$threads" -o "$SORTED" "$BAM"
samtools index "$SORTED"

echo "[se_02] Done."
echo "       Aligned BAM: $SORTED"
