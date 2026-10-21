# 20261021 — 16S reference DB from experiment 20260907, and strain distinguishability

Built as the input for the next experiment: higher-order communities of varying size and
relatedness, read out by nanopore 16S. Analyses `experiments/20260907` only — the run with
two monocultures per strain on different plates plus ~8 pair wells each, so every
sequence can be confirmed by independent wells rather than voted on.

Run `TMPDIR=<scratch> python run_all.py` (env `karl_seq_analysis`, ~9 min; spoa from the
`NGSpeciesID` env, skani from `nasal-strain`, blast from `blast`, barrnap from `NanoSeq`). Fully deterministic — a from-scratch rerun reproduces every output
byte for byte.

## What to use

| file | for |
|---|---|
| `01_reference_db/outputs/reference_db_20260907_haplotypes.fasta` | **read assignment** — one entry per operon haplotype (`<strain>_h<k>`); assign a read to its nearest haplotype, then to the strain |
| `01_reference_db/outputs/reference_db_20260907.fasta` | one majority-base sequence per strain; header carries tier and flags |
| `01_reference_db/outputs/s01_strain_reference.csv` | per strain: tier, evidence, flags, agreement with older DBs |
| `02_distinguishability/outputs/d01_strain_distinguishability.csv` | per strain: identity group, nearest neighbour, how many strains within 0 / 5 / 15 bp |
| `02_distinguishability/outputs/d01_close_pairs.csv` | every pair < 15 bp, with class and operon evidence |
| `02_distinguishability/outputs/d01_identity_groups.csv` | groups of strains with byte-identical 16S |
| `02_distinguishability/outputs/d01_distance_matrix_bp.csv` | full 378×378, haplotype-aware bp |
| `02_distinguishability/outputs/d02_resolution_by_bp_clean.csv` | measured read cross-talk vs separation |
| `02_distinguishability/outputs/d03_max_separable.csv` | largest mutually separable set per threshold |
| **`03_genome_concordance/outputs/g04_design_eligibility.csv`** | **per strain: can it go into a genome-predicted community (`genomic_use` high / medium / exclude, and why), its genome unit** |
| `03_genome_concordance/outputs/g02_strain_genome_status.csv` | per strain: does the genome on file carry the well's 16S, and on what evidence |
| `03_genome_concordance/outputs/g02_genome_clusters.csv` | genomes that are one organism (≥ 99.9% ANI) |
| `03_genome_concordance/outputs/g04_max_separable_genomic.csv` | largest separable community using genome-usable strains only |

## Headline results

**378 / 384 strains have a reference.** Tiers describe independent evidence for what is
in *this run's source plate*:

| tier | n | meaning |
|---|---|---|
| A | 194 | both monocultures, on different plates, give the same sequence |
| B | 102 | one monoculture + ≥2 pair wells |
| C | 31 | no mono but ≥3 pair wells on ≥2 plates (25), or one mono that an older independent run agrees with (6) |
| D | 51 | one mono alone (24), monos disagree with no tiebreak (16), 1–2 pair wells (11) |
| — | 6 | no data: A8 B11 C13 G19 I18 L14 |

22 of the 27 never-grew strains in the DB are tier D. Flags across tiers: `competing` 58
(a second sequence >5 bp away recurs), `mixed_source` 13 (both monos ≥20% second
population), `many_het` 9.

**The consensus is base-exact.** Of 230 strains with two usable monos, 191 (83%) give
identical sequences from the two plates and 199 are within 2 bp; then a gap to 9+ bp,
which is a different organism. The 1–2 bp cases are all at het sites (operons split
~50/50, majority flips between wells). 374 / 378 final sequences match a supporting
observation exactly.

**Independent check.** The older `full_collection` / `contamination_informed` DBs (earlier
runs, other plate copies) agree within 10 bp with 83% of tier A, 85% of B, 87% of C — and
4% of D. The ~15% of A/B that disagree are confirmed by this run's own monos *and* pair
wells, so they describe this source plate; they are likely plate-copy drift.

**Positive control.** The collection sheet lists 9 isolates placed in 2–3 wells each.
11 / 15 of those well pairs come out byte-identical. The 4 misses: MYb71's M24 is 302 bp
from its two siblings (different organism in that well); MYb56's F24 is tier D; MYb57's
H2 and M3 are 17 bp apart.

### Identity

| | |
|---|---|
| pairs with byte-identical 16S | 819, in **72 identity groups covering 254 of 378 strains** |
| largest groups | 27, 17, 12, 10, 7, 7, 7 |
| strains with any partner < 5 bp / < 15 bp | 318 / 347 |
| strains whose nearest neighbour is ≥ 15 bp | 31 |

**Identical 16S is not the same strain**, and operon variation can tell them apart. For each
identical pair the minor-allele frequency at every position is compared (from the reads
behind each reference): 120 pairs differ at some site by ≥ 0.2 → **different genomes**;
291 share the same het pattern; 408 have no het sites (16S has nothing more to say). The
0.2 cut is twice the largest gap seen between known copies of one isolate (0.11).
- **G01** — 27 wells (A10 J10 J6 K14 K16 K18 K20 K22 K24 L12 M2 M6 M8 N2–N22 even, P16 P19
  P22 P24), 119 same-pattern pairs vs 1 different, concentrated in one block of the plate
  (rows K–P, mostly even columns). Most likely one organism in all of them.
- 30 groups (72 strains) are all-same-pattern; 20 groups (109 strains) contain at least one
  different-genome pair; 22 groups (73 strains) are uninformative.

### How close is too close for a nanopore read? Measured, not assumed

For each tier-A strain, reads from one monoculture are raced against a reference from its
other monoculture and against every strain within 40 bp. Clean strains, haplotype-aware
assignment (`d02_resolution_by_bp_clean.csv`):

| separation | reads correct | median cross-talk | p90 | p90, single consensus |
|---|---|---|---|---|
| 0 bp | 11% (87% ties) | 50% | 50% | 50% |
| 1–2 bp | 83% | 5% | 37% | 50% |
| 3–4 bp | 92% | 0.7% | 24% | 23% |
| 5–6 bp | 97.8% | 0 | 3.3% | 10.7% |
| 7–8 bp | 99.4% | 0 | 0 | 0 |
| 9–40 bp | 99.7–99.8% | 0 | ~0.7% | ~0.7% |

From ~7 bp the cross-talk is a flat ~0.7% floor that does **not** improve with distance —
chimeras and stray reads, not resolution. So "<5 bp = closely related" holds; but 15 bp
is conservative — **7–8 bp performs like 15–40 bp**. The earlier 10/15 bp rules came
from 20260630, with references built from a handful of reads and impure monos.

Residual problems above 7 bp are strain-specific, not distance-specific: likely
mixtures (H14/H8, E23, D7), and a few strains (A7, K16, P13, P16) leaking a steady 2–3.5%
to many strains 20–38 bp away — a low-level contaminant in their wells.

### What it means for community size

Largest set with every pair ≥ T bp apart (`d03_max_separable.csv`; exact ILP). *core* =
tier A/B, grew, not `mixed_source`, haplotypes within 10 bp (n = 280); *extended* adds
tier C (n = 307):

| min separation | core | extended |
|---|---|---|
| 15 bp | 85 | 88 |
| 10 bp | 99 | 106 |
| 8 bp | 103 | 110 |
| 7 bp | 105 | 113 |
| 5 bp | 112 | 120 |

## Do the genomes describe what is in the wells? (`03_genome_concordance`)

The genomic tables come from assemblies of the original isolates (panX's 298 GenBank
files). The working plate has drifted. Each well's 20260907 16S (consensus and haplotypes)
was blasted against all 298 assemblies, and the assemblies compared with each other
(skani ANI) and checked for purity (barrnap + SPAdes contig coverage).

**The 298 genomes are 181 distinct organisms.** 162 assemblies fall into 45 clone clusters
(≥ 99.9% ANI; min within-cluster ANI 99.97%, alignment fraction ≥ 94%), the largest of
16. For genome-based ML these are one strain each — `cv_strain` must hold clone-mates out
together, or it leaks.

**Verdict per well** (378 with a reference):

| | ok | wrong | unverified | no genome |
|---|---|---|---|---|
| all | 220 | 61 | 13 | 84 |
| tier A | 128 | 19 | 4 | 43 |

`ok` = the well's 16S is in its own genome (218), or a clone-mate's (2). `wrong` = the
genome carries a clearly different 16S (20), only a contaminant operon matches (3), a
clone-mate's 16S differs (3), or the genome is < 80% ANI to every genome that does carry the
well's 16S (35 — validated: on wells with a known answer this rule flags 15/16 mismatches
and 2/158 matches). **About one in five wells with a genome holds a different organism
from that genome**, and that includes 19 tier-A strains whose two monos agree perfectly —
the monos agree because the source well itself was taken over.

**G01, the 27-well identity group, is one known genome.** 16 wells are the C001 clone
(K14 K16 K18 K20 M2 M6 M8 N14 N16 N18 N2 N20 N22 N6 P22 P24) and their assemblies carry
the G01 16S exactly. The other 9 (A10 J10 J6 L12 N8 N10 N12 P16 P19) hold C001 by 16S but
their genomes are other organisms (A10's is the G02 clone; N8/N10/N12 share another) —
overgrown by C001. K22, K24 have no genome. Overgrowth is concentrated: C001 took 9 wells,
the C002 clone (G02's organism) 5, C034 4.

**9% of assemblies are themselves mixed** (23/294 carry a second organism's 16S). The
tell-tale case: E21 and E22 share an identical full-length 16S but only ~4.8 kb of sequence —
one rRNA operon. In E22's assembly that operon sits at 0.21× bulk coverage, its own 16S at
9.3×: a contaminant in the sequenced culture whose multi-copy operon assembled while the rest
of its genome did not. Matches to a contaminant operon are not counted as `ok`. 16 of 23
mixed assemblies resolve by coverage; 7 stay ambiguous.

**What is left for genome-based community design** (`g04`): 204 strains = **133 distinct
genomes**. `high` (107): no other genome on file carries the same 16S, so the readout pins
down the genome. `medium` (97): consistent, but 16S-identical relatives with different
genomes exist, so a swap would be invisible. Largest separable community:

| min separation | high only | high + medium |
|---|---|---|
| 15 bp | 61 | 77 |
| 8 bp | 69 | 91 |
| 5 bp | 73 | 98 |

## Method, and the traps it avoids

1. **Per-population consensus** (`b01`, `scripts/deep_consensus.py`). spoa over the 50
   highest-qscore reads of each population. Reads are R10.4.1 SUP v5.2, median Q24. On
   P15's two monos, spoa over the top 20/50/100 reads gave the same sequence in both wells;
   medaka moved one by 1 bp (less reproducible, not used); a random 50 including low-Q
   reads drifted 6 bp.
2. **Pair wells, seeded** (`b02`). Greedy read clustering only separates organisms more
   than ~60 bp apart (a Q24 read sits a median ~12 bp from its own consensus), and pairs
   were ≥15 bp apart. So reads go to the nearest candidate sequence of either strain and
   each group gets its own consensus. A group counts only if its consensus is within 15 bp
   of the strain's seeds and closer to them than to the partner's (under competitive
   exclusion a minority's group is mostly the winner's misassigned reads).
3. **One reference per strain** (`b03`). All observations are clustered at ≤2 bp; the
   cluster seen in the most distinct wells wins, monos break ties, and older DBs break
   remaining ties. Het sites count only if they recur in ≥2 independent wells:
   honouring a single well's het sites lets a 50/50 mixed well match either organism at
   0 bp (A17's monos: 21 bp plain, 2 bp het-aware).
4. **Haplotypes** (`b03`). Het sites are linked into operon haplotypes, and a minority
   haplotype can sit much closer to another strain than the consensus does. G15 (5 het
   sites, clean monos) sent 22–31% of its reads to strains 10–11 bp from its consensus;
   with haplotype references that drops to 0–2%.
5. **Distances** (`d01`) ignore 20 bp at each end and take the smaller of the two infix
   directions. Global alignment reported up to 9 bp "difference" between two monos of one
   strain when one consensus kept a few bases of 27F primer tail; taking the larger
   direction put a 19–23 bp floor under every comparison with the older DBs, whose
   sequences still carry primer tails.

Genome-side traps: comparing two partial 16S fragments (or each against a reference) as if
they overlapped invented up to 87 "mixed" assemblies — genes are compared only where they
actually overlap; the concatenated BLAST FASTA must not sit where skani globs genomes (it
matched everything); a 400–800 bp conserved stretch matches across distant taxa, so a genome
"carries" a 16S only over ≥ 1000 bp.

Bugs found and fixed on the way: `het_sites` recorded a "minor" allele equal to the final
consensus base at ~50/50 sites (G15 position 431 as A/A), which dropped those sites from
phasing and from the IUPAC reference; a "top minor base only" profile called MYb181's two
wells different genomes because position 180 is tri-allelic.

## Not done / open

- Taxonomy is not attached. The sheet's `phylogeny` column belongs to the original-isolate
  lineage, and only 219 / 384 wells have an identity agreed across that lineage and the
  working plates (`merge_consensus_sequences/external_reference_crosscheck`). Classifying these
  sequences directly (SILVA/GTDB) would give a usable relatedness axis.
- The genomic ML built on 20260907 (`analyses/20260921/07_genomic_ml`) joined by well name,
  so it trained on a wrong genome for ~1 in 5 joined wells, and its `cv_strain` folds let
  clone-mates straddle train and test. Not yet re-run with `genomic_use` / `genomic_unit`.
  20260630 (the plate as it was in July) has not been checked this way; `g01`/`g02` would
  do it given that run's references.
- The DB describes the **20260907 source plate**. A fresh copy from the freezer master may
  differ for some wells; include monocultures to re-anchor.
