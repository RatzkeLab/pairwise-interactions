# Analyses using experiment 20260721

One row per analysis run that touched this experiment. Read this instead of opening the
analysis folders; open them only for method detail or for code worth reusing.

> **Read [known_issues.md](known_issues.md) first.** This experiment's well labels do not
> identify the organisms in the wells. Every result below is conditioned on that.

| analysis | date | what it asked | what it found |
|---|---|---|---|
| [`analyses/20260804`](../../analyses/20260804) | 2026-08 → 09 | Reference DB, read QC, relative abundance, ML join test | Built an 81-strain consensus and the interaction pipeline. The genomic join test **failed**, which is what first exposed the label problem. |
| [`analyses/20260827`](../../analyses/20260827) | 2026-08 → 09 | Is this the same source plate as 20260630? Can the labels be recovered? | OD says **yes, same stock** — so the mix-up is downstream of the freezer. About half the labels are recoverable by 16S, validated against interaction labels rather than circularly. |
| [`analyses/20260921`](../../analyses/20260921) | 2026-09 | Used as the **negative control** for 20260907's identity check | 0 supported of 188 decidable comparisons, across all three reference databases. |

## This experiment's standing

**Usable for method development and as a negative control; not usable for genomic ML.**
The sequencing is fine — the reads are real and the consensus sequences are
self-consistent. What fails is the mapping from well label to organism.

The label recovery in `analyses/20260827` is the route back to using this data, and it
would benefit from a re-demultiplex (see [demux_runs.md](demux_runs.md)).
