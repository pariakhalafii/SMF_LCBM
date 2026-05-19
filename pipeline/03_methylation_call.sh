#!/usr/bin/env bash
# 03_methylation_call.sh -- bismark_methylation_extractor with NOMe-seq / SMF flags.
#
# Produces:
#   * bedGraph + coverage files (CpG and non-CpG)
#   * a CX_report.txt.gz with per-cytosine, context-resolved methylation calls
#
# This is a *bulk* methylation sanity check ("X% of all Cs are methylated").
# For real per-read SMF analysis use scripts/extract_per_read_methylation.py.
#
# Usage: bash pipeline/03_methylation_call.sh config/config.example.yaml
#
# Requires: bismark, samtools, python3.

set -euo pipefail

CONFIG="${1:-config/config.example.yaml}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

eval "$(python3 "$HERE/_load_config.py" "$CONFIG")"

: "${sample_id:?missing sample_id}"
: "${genome_dir:?missing genome_dir}"
: "${out_dir:?missing out_dir}"
: "${threads:=4}"
: "${extract_options:=--CX --cytosine_report --comprehensive}"

ALIGN_DIR="$out_dir/align"
METH_DIR="$out_dir/methylation"
mkdir -p "$METH_DIR"

# Pick the most-processed BAM that exists, in priority order.
# Paired-end Bismark output has a "_pe" suffix; single-end output (e.g. merged
# FLASH reads) does not.
BAM=""
MODE=""
for pattern in "*_pe.sorted.bam:paired" \
               "*_pe.bam:paired" \
               "*minlen*.sorted.bam:single" \
               "*deduplicated.bam:auto" \
               "*_bismark_bt2.bam:single"; do
    glob="${pattern%%:*}"
    m="${pattern##*:}"
    cand=$(ls "$ALIGN_DIR"/$glob 2>/dev/null | head -n1 || true)
    if [[ -n "$cand" && -f "$cand" ]]; then
        BAM="$cand"
        MODE="$m"
        break
    fi
done

if [[ -z "$BAM" ]]; then
    echo "ERROR: no Bismark BAM found in $ALIGN_DIR" >&2
    echo "       Run 02_align.sh first." >&2
    exit 1
fi

# Auto-detect paired-end vs single-end if needed (deduplicated.bam keeps suffix
# from upstream, so we trust the original name).
if [[ "$MODE" == "auto" ]]; then
    if [[ "$BAM" == *"_pe"* ]]; then
        MODE="paired"
    else
        MODE="single"
    fi
fi

if [[ "$MODE" == "paired" ]]; then
    PE_FLAG="--paired-end"
    echo "[03] Detected paired-end BAM"
else
    PE_FLAG="--single-end"
    echo "[03] Detected single-end BAM"
fi

echo "[03] Using BAM: $BAM"
echo "[03] Running bismark_methylation_extractor ($PE_FLAG)"
bismark_methylation_extractor \
    $PE_FLAG \
    --bedGraph \
    --gzip \
    --multicore "$threads" \
    --genome_folder "$genome_dir" \
    $extract_options \
    -o "$METH_DIR" \
    "$BAM"

echo "[03] Done. CX_report and bedGraph files in $METH_DIR"
echo "[03] For per-read SMF analysis: scripts/extract_per_read_methylation.py"
