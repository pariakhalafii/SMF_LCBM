# What the code does, step by step:
Step 1 — Clean the raw reads. It runs FastQC (a quality report) and Trim Galore (chops off the sequencing adapters and any low-quality bases at the ends). Output: cleaner FASTQs in results/<sample>/trimmed/.


Step 1b — Glue R1 and R2 together with FLASH. Your two reads overlap in the middle, so FLASH stitches each pair into one long read covering the whole molecule. You set "min overlap 30, max overlap 100." Output: a single FASTQ of merged reads in results/<sample>/merged/, plus a histogram of how long the merged reads are.



Step 2 — Align with Bismark. It compares each merged read to your reference and figures out where it came from, while also accounting for bisulfite/EM-seq C→T conversion. It saves the reads that didn't align so you can sanity-check them. Output: a BAM file (the alignment) in results/<sample>/align/.




Step 3 — Sort and index the BAM. This is just bookkeeping so IGV (and the rest of the pipeline) can jump to any region quickly. Output: a .sorted.bam and a .bai index file next to it.




Step 4 — You eyeball it in IGV. Not code — you open the sorted BAM in IGV and check that reads land where you expect, look full-length, show C→T changes, etc. Skip-able if you trust your sample.



Step 5 — Drop short reads. Anything shorter than 500 bp gets thrown out, because short fragments leave gaps in your heatmap and bias the average. Output: a .minlen500.sorted.bam.




Step 6 — Bulk methylation sanity check (Bismark again). Tells you "X% of all your Cs are methylated" overall. You're not analyzing this — you're just checking that the experiment worked. Output: a CX_report.txt.gz and bedGraphs in results/<sample>/methylation/.




Step 7 — Custom per-read methylation caller. This is the real analysis. For every read, it walks along and asks at each informative C: "is this a GpC (accessibility) or a CpG (endogenous)? Is it methylated or not?" It throws out the ambiguous GCG sites. Output: a TSV table where every row is one C on one read
Plus an optional "molecule pattern" file with one row per read showing the whole pattern as a string of M/U/. characters.





Step 8 — Plots. Five outputs in examples/ (or wherever you point it):

Single-molecule heatmap — one row per read, colored cells for methylated/unmethylated/no-call. Sorted so the most-protected molecules cluster together. This is the SMF picture.
Average accessibility curve — a smooth line showing the average GpC methylation across the window. Dips where things are protected (TF, nucleosome).
State composition bar — what fraction of your molecules are accessible vs. TF-bound vs. nucleosome.
Coverage per position — how many reads cover each C. Red bars = positions you should not trust.
Fragment length distribution — histogram of merged-read lengths, with your length cutoff drawn on it.

Plus a CSV summary table and a confusion matrix (true vs. predicted state, useful when validating).