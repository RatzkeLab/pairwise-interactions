"""Constructed metabolic pair-features: resource competition and cross-feeding potential.

Everything tried so far fed the model *generic* per-strain feature vectors (sum, difference) and
let it work out the interaction. These features instead encode specific ecological mechanisms
directly as pair-level quantities.

**BiGG cannot supply these, despite the name.** `BiGG_and_strains_table.csv` has 9442 columns of
the form `<model_id>.<gene_locus>` -- 73 published genome-scale models (mostly E. coli, Shigella,
Klebsiella variants) x their gene loci. They are gene-orthology calls, NOT reactions, metabolites,
exchange fluxes or stoichiometry, so no uptake overlap, no byproduct exchange and no joint FBA
can be computed from that file. Genuine flux-based features would require building draft GEMs
from the assemblies (CarveMe / gapseq / ModelSEED) and running something like SMETANA -- a real
bioinformatics step, not a transformation of the tables we have. That also explains why BiGG was
among the worst feature tables in the sweep: for a collection of Pseudomonas/Sphingomonas it is
an E. coli-centric orthology table, and a redundant one (the same locus appears under up to 4
models).

What the available tables CAN support:

  **Resource competition -- CAZy.** Carbohydrate-active enzyme families split by prefix into
  degradative (GH glycoside hydrolases, PL lyases, CE esterases, CBM binding modules) versus
  biosynthetic (GT glycosyltransferases). Two strains sharing degradative families can attack
  the same carbohydrates, so overlap there is a direct competition proxy. The split is by
  prefix, so it needs no external annotation.

  **Metabolic complementarity -- KEGG Modules.** Pathway-level presence/absence: shared modules
  (redundancy), modules unique to each partner (complementarity), and the asymmetry between
  them.

  **Cross-feeding potential -- KEGG biosynthesis modules.** If A lacks a biosynthesis module
  that B has, B can potentially feed A. Directional, so it gives a genuinely antisymmetric
  feature for the who-wins target and a summed one for yield.

    `BIOSYNTHESIS_MODULES` was VERIFIED against https://rest.kegg.jp/list/module on 2026-08-28.
    The first draft, written from memory, contained five wrong entries (a catabolic shunt, two
    C1-interconversion modules, and two eukaryote-only routes); those are removed and listed in
    the comment beside the constant. Results are still reported with and without the
    cross-feeding block, because that block carries a large share of the signal and its
    contribution should stay visible rather than being folded into a single number.
"""

import numpy as np
import pandas as pd

# Amino-acid and cofactor/vitamin biosynthesis modules.
# VERIFIED 2026-08-28 against the KEGG REST API (https://rest.kegg.jp/list/module, 573 modules):
# every ID below is present in KEGG and its name contains "biosynthesis", except M00022
# (shikimate pathway), kept deliberately -- it produces chorismate, the obligate precursor of
# all three aromatic amino acids, so losing it is an aromatic-AA auxotrophy.
# Removed after verification, having been wrong in the first draft:
#   M00027  GABA shunt              -- catabolic/interconversion, not biosynthesis
#   M00140  C1-unit interconversion -- interconversion, not biosynthesis
#   M00141  C1-unit interconversion (eukaryotes) -- also not a bacterial pathway
#   M00128  Ubiquinone biosynthesis, EUKARYOTES  -- prokaryote route is M00117
#   M00868  Heme biosynthesis, ANIMALS AND FUNGI -- bacterial route is M00121
# M00846/M00880 are genuine biosynthesis modules but absent from this collection's table.
BIOSYNTHESIS_MODULES = [
    # amino acids
    "M00015", "M00016", "M00017", "M00018", "M00019", "M00020", "M00021", "M00022",
    "M00023", "M00024", "M00025", "M00026", "M00028", "M00033", "M00338",
    "M00432", "M00525", "M00526", "M00527", "M00535", "M00570", "M00844", "M00845",
    # cofactors / vitamins (prokaryotic routes only)
    "M00115", "M00116", "M00117", "M00119", "M00120", "M00121", "M00122", "M00123",
    "M00124", "M00125", "M00126", "M00127", "M00577", "M00840", "M00841",
    "M00842", "M00843", "M00846", "M00880",
]

DEGRADATIVE_CAZY_PREFIXES = ("GH", "PL", "CE", "CBM")


def _presence(df):
    return (df > 0).astype(int)


def load_table(path, strains, w2s):
    """Table indexed by well label, restricted to the strains actually modelled."""
    t = pd.read_csv(path, index_col=0)
    t.index = t.index.astype(str)
    genomes = [w2s[s] for s in strains]
    keep = [g for g in genomes if g in t.index]
    sub = t.loc[keep]
    sub.index = [s for s, g in zip(strains, genomes) if g in t.index]
    return _presence(sub)


def build(pairs, cazy=None, modules=None, ko=None, include_crossfeeding=True):
    """Pair-level constructed features.

    Returns (DataFrame indexed like `pairs`, symmetric_cols, antisymmetric_cols). Symmetric
    columns are unchanged when the two strains swap; antisymmetric ones negate. The yield target
    should use the symmetric block only -- an antisymmetric feature cannot inform a quantity
    that does not change when you relabel the partners.
    """
    A, B = pairs.strain_a.values, pairs.strain_b.values
    out, sym, anti = {}, [], []

    def add_pairwise(mat, tag):
        if mat is None or mat.shape[1] == 0:
            return
        a = mat.reindex(A).values
        b = mat.reindex(B).values
        ok = ~(np.isnan(a).any(1) | np.isnan(b).any(1))
        a = np.nan_to_num(a)
        b = np.nan_to_num(b)
        shared = (a * b).sum(1)
        only_a = (a * (1 - b)).sum(1)
        only_b = (b * (1 - a)).sum(1)
        union = np.maximum(((a + b) > 0).sum(1), 1)
        out[f"{tag}_shared"] = shared                       # redundancy / competition
        out[f"{tag}_jaccard"] = shared / union              # niche overlap
        out[f"{tag}_complementary"] = only_a + only_b       # combined unique capability
        out[f"{tag}_total"] = union
        sym.extend([f"{tag}_shared", f"{tag}_jaccard", f"{tag}_complementary", f"{tag}_total"])
        out[f"{tag}_advantage_a_minus_b"] = only_a - only_b  # directional capability asymmetry
        anti.append(f"{tag}_advantage_a_minus_b")
        out[f"{tag}_valid"] = ok.astype(float)
        sym.append(f"{tag}_valid")

    if cazy is not None:
        deg = cazy[[c for c in cazy.columns if c.startswith(DEGRADATIVE_CAZY_PREFIXES)]]
        bio = cazy[[c for c in cazy.columns if c.startswith("GT")]]
        add_pairwise(deg, "cazy_degradative")   # resource competition
        add_pairwise(bio, "cazy_biosynthetic")
    if modules is not None:
        add_pairwise(modules, "kegg_module")
    if ko is not None:
        add_pairwise(ko, "ko")

    if include_crossfeeding and modules is not None:
        bio_cols = [c for c in modules.columns if c in set(BIOSYNTHESIS_MODULES)]
        if bio_cols:
            m = modules[bio_cols]
            a = np.nan_to_num(m.reindex(A).values)
            b = np.nan_to_num(m.reindex(B).values)
            # B has the biosynthesis module, A lacks it -> B can potentially feed A
            b_feeds_a = (b * (1 - a)).sum(1)
            a_feeds_b = (a * (1 - b)).sum(1)
            out["crossfeed_mutual"] = b_feeds_a + a_feeds_b
            out["crossfeed_shared_biosynth"] = (a * b).sum(1)
            sym.extend(["crossfeed_mutual", "crossfeed_shared_biosynth"])
            out["crossfeed_asymmetry_a_self_sufficiency"] = a_feeds_b - b_feeds_a
            anti.append("crossfeed_asymmetry_a_self_sufficiency")
            out["n_biosynthesis_modules_used"] = np.full(len(A), len(bio_cols), float)

    F = pd.DataFrame(out, index=pairs.index)
    sym = [c for c in sym if c in F.columns and F[c].std() > 0]
    anti = [c for c in anti if c in F.columns and F[c].std() > 0]
    return F, sym, anti
