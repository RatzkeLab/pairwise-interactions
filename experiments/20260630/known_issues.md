# Known issues — experiment 20260630

Defects in the data or in the methods applied to it, and what they affect.

---

## <a name="minibar-barcode-file"></a>Malformed barcode file — FIXED 2026-09-24 (re-demultiplexed)

**Found** 2026-09-21 while working on 20260907. **Affects** `20260710_demultiplex` and
every analysis built on it. **Status** fixed: re-demultiplexed as
`20260710_demultiplex_increased_tolerance` (4.18× the reads, false-assignment rate 0.50% →
0.61%), and re-analysed side by side with the old demux in
`analyses/20260710_increased_tolerance`. The genomic-ML headline did not move (see
`analysis_index.md`). `analyses/20260710` itself still reflects the old demux.

minibar reads column 2 of the barcode TSV as the forward index. `minibar_primers_20260630.tsv`
puts the whole 59 bp construct there — 15 bp adapter + 24 bp barcode + 20 bp 16S primer —
so minibar matches a 59-mer within the barcode edit distance instead of the 24-mer
barcode. Same tolerance over 2.5× the length, and most reads fail into `unk.fastq`.

Measured on 20260907, which has the identical file layout: assignment went from 13.5% to
**64.1%** of filtered reads after correcting the file and setting `-e 4`. This experiment
assigns 14.9%, so expect a comparable gain.

**Nothing already assigned is wrong** — the assigned reads were assigned correctly. There
are simply ~4× more of them available. Fix and calibration:
`analyses/20260921/02_demux_qc/{fix_barcode_file.py,calibrate_minibar.py}`.

---

## Single-linkage clustering chains at depth — FIXED in shared code

**Found** 2026-09-22. **Affects** `analyses/20260710/04_qc` (`check_reads.ipynb`) and any
per-well strain-count result from it. **Status** method fixed in
`scripts/well_composition.py`; the old outputs were not regenerated.

Single-linkage connected components merge the two strains of a pair well once there are
enough error-rich reads to bridge the same/different distance modes. It gets *worse* with
more data — on 200 pair wells, wells resolving both strains went 27.0% → 10.0% between 120
and 500 reads, while greedy seed clustering held at 47.0% → 45.0%.

So **two-strain counts from that notebook are understated.** The replacement decides
membership against a cluster seed rather than a neighbour and cannot chain.

---

## Genome-derived 16S is a weak reference — INHERENT

**Affects** any identity check against `rDNA_16S_db_all_strains.fasta`.

Half the entries are assembly fragments, as short as 396 bp against a ~1,450 bp consensus,
so global alignment scores them as different organisms regardless of true identity. Use
overlap-aware distance (`scripts/strain_identity.flexible_distance`): doing so moved this
experiment's label support against that reference from 11.5% to **34.1%**.

Separately, 16S simply cannot separate much of this collection — only 16.2% of the
collection reference's own sequences are distinguishable at 10 bp. A wrong-named best hit
is usually meaningless and must be scored *ambiguous*, not *contradicted*.

---

## Derived data is large and untracked

`analyses/20260710` carries ~39,000 output files. Outputs are gitignored by design, so
they exist only on this machine; regenerating them requires re-running the analysis.
