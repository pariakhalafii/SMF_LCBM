#!/usr/bin/env bash
# 02b_filter.sh -- Length-filter the sorted BAM to drop short fragments.
#
# Reads shorter than your amplicon (e.g. 3N601 constructs) bias the heatmap by
# leaving NaN columns and inflating "accessible" calls in regions covered only by
# truncated molecules. We drop anything below `min_read_length` (default 500 bp).
#
# Usage: bash pipeline/02b_filter.sh config/config.example.yaml
#
# Requires: samtools, python3.

set -euo pipefail

CONFIG="${1:-config/config.example.yaml}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(python3 "$HERE/_load_config.py" "$CONFIG")"

: "${sample_id:?missing sample_id}"
: "${out_dir:?missing out_dir}"
: "${threads:=4}"
: "${min_read_length:=500}"

ALIGN_DIR="$out_dir/align"
INPUT_BAM=$(ls "$ALIGN_DIR"/*_bismark_bt2*.sorted.bam | head -n1)
FILT_BAM="${INPUT_BAM%.sorted.bam}.minlen${min_read_length}.sorted.bam"

echo "[02b] Filtering $INPUT_BAM for read length >= $min_read_length"
# samtools view -e supports an inline expression; "length(seq)" gives the read length
# without needing to compute CIGAR-based reference length.
samtools view -@ "$threads" -h -b -e "length(seq)>=${min_read_length}" "$INPUT_BAM" -o "$FILT_BAM"
samtools index "$FILT_BAM"

# Quick before/after counts so the user can see how many reads survived.
N_IN=$(samtools view -c "$INPUT_BAM")
N_OUT=$(samtools view -c "$FILT_BAM")
PCT=$(python3 -c "print(f'{100*${N_OUT}/${N_IN}:.1f}' if ${N_IN} else 'NA')")
echo "[02b] Kept $N_OUT / $N_IN reads (${PCT}%)"
echo "[02b] Filtered BAM: $FILT_BAM"
