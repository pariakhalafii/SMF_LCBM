# Single Molecule Footprinting (SMF) — bisulfite / EM-seq, Illumina

A Python toolkit + shell pipeline for **dual-enzyme SMF** (GpC + CpG methyltransferase
treatment, bisulfite or EM-seq conversion, paired-end Illumina sequencing) in the style
of Krebs / Sönmezer et al.

The pipeline goes from raw FASTQs all the way to single-molecule plots and runs both on
real data and on a built-in synthetic example.

---

## Protocol overview

The shell pipeline mirrors the eight-step protocol below:

| # | Step | Tool | Script |
|---|------|------|--------|
| 1 | QC + adapter trim of paired-end FASTQs | FastQC, Trim Galore | `pipeline/01_preprocess.sh` |
| 1b | **Merge overlapping R1/R2** (paired reads overlap in the middle) | **FLASH** (min overlap 30, max 100) | `pipeline/01b_merge.sh` |
| 2 | Bismark genome prep + alignment of merged reads (single-end, `--non_directional`, save `--un` / `--ambiguous`) | Bismark | `pipeline/02_align.sh` |
| 3 | Sort + index BAM for IGV (built into step 2) | samtools | `pipeline/02_align.sh` |
| 4 | Visual inspection in IGV | IGV | manual |
| 5 | **Length filter** (≥ 500 bp) on the sorted BAM to drop short fragments | samtools | `pipeline/02b_filter.sh` |
| 6 | Bulk methylation as sanity check | `bismark_methylation_extractor` | `pipeline/03_methylation_call.sh` |
| 7 | **Custom per-read methylation caller** (read_id × position table) | smf + pyfaidx | `scripts/extract_per_read_methylation.py` |
| 8 | Plots: heatmap, average accessibility, coverage per position, fragment length, state composition | matplotlib | `smf.viz` + `scripts/run_example.py` |

### Why FLASH and not bwa pemerge

bwa pemerge is unreliable on bisulfite/EM-seq-converted reads (lots of C→T mismatches in
the overlap break its consensus). FLASH is much more forgiving and exposes the
`--min-overlap` / `--max-overlap` knobs you need for ~50 bp overlaps. Defaults in
`config/config.example.yaml` are min 30, max 100.

### Why Bismark and not bwa-meth / BSBolt

You can align bisulfite/EM-seq reads with `bwa-meth` or `BSBolt` and they're noticeably
faster, but for SMF specifically Bismark is the right choice because it writes a
**per-base methylation call string into the BAM `XM` tag** for every read. That's what
the custom per-read caller in step 7 consumes — without it you'd have to round-trip
through `MethylDackel` and reassemble per-read state from per-cytosine bedGraphs. Stick
with Bismark.

EM-seq reads align with Bismark exactly the same way as bisulfite reads (`--non_directional`).
Nothing else in the downstream analysis changes.

### Why a custom per-read caller

`bismark_methylation_extractor` produces *bulk* per-cytosine averages — useful as a
sanity check ("23% methylation across all Cs") but not what we need. SMF is a
single-molecule readout, so we need a binary call at every informative cytosine on
every read, partitioned correctly between **GpC** (M.CviPI accessibility readout) and
**CpG** (endogenous methylation), with **GCG** ambiguous sites excluded. Bismark calls
CpG/CHG/CHH, but it doesn't natively split GpC out from CHH/CHG. The custom caller in
`scripts/extract_per_read_methylation.py` resolves context from the reference and
emits the per-read table.

### Dedup: amplicon vs. whole-genome

For amplicon SMF (e.g. 3N601 nucleosome-positioning constructs) every read is at the
same locus by design, so positional `deduplicate_bismark` will collapse legitimately
distinct molecules. **Dedup is OFF by default** in `config/config.example.yaml`. Set
`deduplicate: yes` only for whole-genome SMF.

---

## Layout

```
SMF/
├── README.md
├── requirements.txt
├── config/
│   └── config.example.yaml         # Sample config for the shell pipeline
├── pipeline/
│   ├── 01_preprocess.sh            # FastQC + Trim Galore
│   ├── 01b_merge.sh                # FLASH merge of overlapping R1/R2
│   ├── 02_align.sh                 # Bismark alignment (merged SE, --non_directional)
│   ├── 02b_filter.sh               # samtools length filter (>= 500 bp by default)
│   └── 03_methylation_call.sh      # bismark_methylation_extractor (sanity check)
├── smf/
│   ├── __init__.py
│   ├── context.py                  # GpC / CpG / GCG context detection
│   ├── io.py                       # CX reports + per-read BAM XM parsing
│   ├── per_read.py                 # Per-read methylation matrices
│   ├── classify.py                 # Single-molecule footprint classification
│   ├── stats.py                    # Co-occupancy & average accessibility
│   └── viz.py                      # Heatmap, average, coverage, fragment-length plots
├── scripts/
│   ├── extract_per_read_methylation.py   # Custom per-read caller (step 7)
│   ├── make_synthetic_data.py            # Synthetic per-read matrix
│   └── run_example.py                    # End-to-end downstream demo
├── tests/
│   ├── test_smf.py
│   └── run_tests.py                # stdlib-only runner; pytest also works
└── examples/                        # Plots dropped here when run_example.py runs
```

## Install

The Python package only needs `numpy`, `pandas`, `matplotlib`. `scipy` is reserved for
future use; `pysam` and `pyfaidx` are needed only for the BAM/FASTA-based caller in
step 7. `pytest` is optional (a stdlib runner is included).

```bash
pip install -r requirements.txt
```

## Running the upstream pipeline (real data)

Edit `config/config.example.yaml` for your sample, then:

```bash
bash pipeline/01_preprocess.sh        config/config.example.yaml   # FastQC + Trim Galore
bash pipeline/01b_merge.sh            config/config.example.yaml   # FLASH (min 30 / max 100)
bash pipeline/02_align.sh             config/config.example.yaml   # Bismark + sort + index
bash pipeline/02b_filter.sh           config/config.example.yaml   # length >= 500 bp
bash pipeline/03_methylation_call.sh  config/config.example.yaml   # bulk methylation sanity check
```

Then run the per-read caller and produce per-read calls + molecule-pattern strings:

```bash
python scripts/extract_per_read_methylation.py \
    --bam   results/demo/align/*.minlen500.sorted.bam \
    --fasta refs/genome/genome.fa \
    --region chr1:10000-12000 \
    --contexts GCH,HCG \
    --out                results/demo/methylation/per_read_calls.tsv \
    --molecule-pattern   results/demo/methylation/molecule_patterns.tsv
```

`per_read_calls.tsv` has one row per (read, informative cytosine):

```
read_id  chrom  pos  strand  ref_base  context  call          methylated
read_1   chr1   10042  +       C        GCH      methylated     1
read_1   chr1   10058  +       C        HCG      unmethylated   0
...
```

`molecule_patterns.tsv` has one row per read with a string pattern (`M`/`U`/`.`) across
all informative cytosines.

## Sanity-check checklist (step 4, in IGV)

Open the sorted BAM in IGV and confirm:

* Reads align to the expected reference.
* You see C→T changes (bisulfite/EM-seq conversion is happening).
* Reads are full length (≈ amplicon size).
* No unexpected pile-ups of short fragments.
* FP sample looks visibly different from the NoFP sample.

## Running the downstream demo (no real data needed)

```bash
python scripts/run_example.py
```

Synthesises 600 reads across a 600 bp window centred on a "TF motif" with three
populations (accessible / TF-bound / nucleosome) and produces, in `examples/`:

* `single_molecule_heatmap.png` — sorted single-molecule heatmap
* `average_accessibility.png` — smoothed bulk readout (drops at the motif)
* `state_composition.png` — % accessible / TF-bound / nucleosome
* `coverage_per_position.png` — reads per informative cytosine, with low-coverage
  positions highlighted in red
* `fragment_length.png` — fragment-length distribution with median + length filter
* `per_read_summary.csv`, `confusion_matrix.csv` — numeric tables

## Tests

```bash
# pytest if installed
pytest -q

# or the stdlib runner
python tests/run_tests.py
```

## Biology, briefly

In dual-enzyme SMF:

* **M.CviPI** methylates accessible **GpC** sites — GpC methylation reports
  *chromatin accessibility* on each individual molecule.
* Endogenous **CpG** methylation is preserved through bisulfite/EM-seq conversion and
  reports the cell's native methylome.
* **GCG** sites are ambiguous (the C is part of both contexts) and must be excluded.
* Each merged read covers tens of GpC sites, giving a binary accessibility vector per
  molecule across the window.

The standard analysis classifies each molecule at a TF binding site into one of:

| State              | Footprint signature                                |
|--------------------|----------------------------------------------------|
| Fully accessible   | High GpC methylation across the whole window       |
| TF-bound           | Small (~30 bp) protected patch over the motif      |
| Nucleosome         | Large (~150 bp) protected patch covering the motif |
