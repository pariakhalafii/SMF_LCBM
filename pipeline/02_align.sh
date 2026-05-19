#!/usr/bin/env bash
# 02_align.sh -- Bismark paired-end alignment for SMF.
#
# Workflow:
#   * Bismark genome preparation (one-time, cached)
#   * Bismark PAIRED-END alignment of trimmed R1/R2 reads (--non_directional)
#       saving --un and --ambiguous so you can inspect what failed to align
#   * Optional deduplication (off by default for amplicon SMF -- see config)
#   * samtools sort + index so the BAM is IGV-ready
#
# Usage: bash pipeline/02_align.sh config/config.example.yaml
#
# Requires: bismark, bowtie2 (or hisat2), samtools, python3.
#
# NOTE: This is the paired-end version of the pipeline. Use it when R1/R2 don't
#       overlap much (typical real-world libraries). For heavily overlapping
#       reads on a short amplicon, you can FLASH-merge first with 01b_merge.sh
#       and align the merged single-end FASTQ -- but that path is rarely needed.

set -euo pipefail

CONFIG="${1:-config/config.example.yaml}"
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

# Trim Galore appends _val_1 / _val_2 to the paired output.
R1_TRIM=$(ls "$TRIM_DIR"/*_val_1.fq.gz 2>/dev/null | head -n1)
R2_TRIM=$(ls "$TRIM_DIR"/*_val_2.fq.gz 2>/dev/null | head -n1)

if [[ -z "$R1_TRIM" || -z "$R2_TRIM" ]]; then
    echo "ERROR: trimmed paired-end FASTQs not found in $TRIM_DIR" >&2
    echo "       Run pipeline/01_preprocess.sh first." >&2
    exit 1
fi

if [[ ! -f "$genome_dir/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2" ]]; then
    echo "[02] Bismark genome preparation (one-time)"
    bismark_genome_preparation --bowtie2 "$genome_dir"
fi

echo "[02] Bismark paired-end alignment"
echo "       R1:      $R1_TRIM"
echo "       R2:      $R2_TRIM"
echo "       options: $bismark_options"
bismark \
    --genome "$genome_dir" \
    -1 "$R1_TRIM" -2 "$R2_TRIM" \
    -p "$threads" \
    $bismark_options \
    -o "$ALIGN_DIR"

# Bismark paired-end output is *_val_1_bismark_bt2_pe.bam
BAM=$(ls "$ALIGN_DIR"/*_pe.bam | head -n1)

# Use `tr` to lowercase rather than ${var,,}, which requires bash 4+ (macOS ships 3.2).
DEDUP_LC=$(printf '%s' "$deduplicate" | tr '[:upper:]' '[:lower:]')
if [[ "$DEDUP_LC" == "yes" || "$DEDUP_LC" == "true" ]]; then
    echo "[02] Deduplicating $BAM (whole-genome SMF mode)"
    deduplicate_bismark --paired --bam --output_dir "$ALIGN_DIR" "$BAM"
    BAM=$(ls "$ALIGN_DIR"/*deduplicated.bam | head -n1)
else
    echo "[02] Skipping deduplication (amplicon SMF mode -- set deduplicate: yes for WG)"
fi

SORTED="${BAM%.bam}.sorted.bam"
echo "[02] Sorting & indexing $SORTED for IGV / random access"
samtools sort -@ "$threads" -o "$SORTED" "$BAM"
samtools index "$SORTED"

echo "[02] Done."
echo "       Aligned BAM:        $SORTED"
echo "       Unaligned reads:    $ALIGN_DIR/*_unmapped_reads.fq.gz"
echo "       Ambiguous reads:    $ALIGN_DIR/*_ambiguous_reads.fq.gz"
echo
echo "[02] Sanity-check in IGV:"
echo "       * Are reads aligned to the expected reference?"
echo "       * Do you see C/T changes (bisulfite/EM-seq conversion)?"
echo "       * Are mate pairs both mapped?"
echo "       * Mapping efficiency in the *_PE_report.txt should be >70% for clean libraries."
