"""Wells -> forward x reverse barcode combinations -> minibar sample sheet.

The primer source plate is a 384 checkerboard: 96 forward primers in rows A,C,E,..,
columns 1,3,5,..; 96 reverse primers in rows B,D,..., columns 2,4,.... Every forward x
reverse pair with unequal index is a valid combination (9,120). Each ordered oligo in
PrimerPlateSpecs is 15 bp ligation adapter + 24 bp barcode + 20 bp 16S primer; its rows
are the 96 reverse primers (in reverse-source-well order), then the 96 forward.

THE BUG THIS MODULE EXISTS TO PREVENT: 20260630, 20260721 and 20260907 all wrote the whole
59 bp construct into minibar's index column, which costs ~5x the reads. Here the index
column is the 24 bp barcode, cut exactly as `fix_barcode_file.py` cut it (after asserting
the construct ends in the primer), and `validate.check_minibar` refuses anything else.
"""
import numpy as np
import pandas as pd

from . import plates

ADAPTER_LEN = 15
BARCODE_LEN = 24
FWD_PRIMER = "AGRGTTYGATYMTGGCTCAG"
REV_PRIMER = "CGGYTACCTTGTTACGACTT"


def primer_source_wells():
    grid = np.array([[f"{r}{c}" for c in range(1, plates.N_COLS + 1)] for r in plates.ROWS])
    return grid[::2, ::2].flatten(), grid[1::2, 1::2].flatten()


def barcode_table(primer_specs_csv):
    """side, idx, source_well, barcode (24 bp), full oligo."""
    fwd_wells, rev_wells = primer_source_wells()
    s = pd.read_csv(primer_specs_csv, usecols=["Well Position", "Sequence Name", "Sequence"])
    s["oligo"] = s.Sequence.str.replace(" ", "", regex=False)
    if len(s) != 192:
        raise ValueError(f"expected 192 primers, found {len(s)}")
    rows = []
    for i, r in enumerate(s.itertuples()):
        side, idx = ("rev", i) if i < 96 else ("fwd", i - 96)
        primer = REV_PRIMER if side == "rev" else FWD_PRIMER
        if not r.oligo.endswith(primer):
            raise ValueError(f"{r._2}: oligo does not end in the {side} primer")
        bc = r.oligo[ADAPTER_LEN:-len(primer)]
        rows.append(dict(side=side, idx=idx, name=r._2,
                         source_well=(rev_wells if side == "rev" else fwd_wells)[idx],
                         barcode=bc, oligo=r.oligo))
    t = pd.DataFrame(rows)
    if (t.barcode.str.len() != BARCODE_LEN).any():
        raise ValueError("barcode length is not 24 bp after stripping adapter and primer")
    return t


def failed_barcodes(fwd_csv, rev_csv, q_max=0.05):
    """Barcodes whose wells dropped out in 20260907 far beyond chance: the oligo, not the
    sample (their wells had cells by OD, and the same barcodes worked in earlier runs).
    Same criterion as analyses/20260921/02_demux_qc/barcode_forensics.py."""
    out = {}
    for side, path in (("fwd", fwd_csv), ("rev", rev_csv)):
        f = pd.read_csv(path)
        base = f.k.sum() / f.n.sum()
        out[side] = set(f[(f.q < q_max) & (f.rate > base)][f"{side}_idx"].astype(int))
    return out


def assign(wells, rng, exclude=None):
    """wells: list of (barcode_plate, well). Each gets a unique (fwd, rev) combination."""
    exclude = exclude or {"fwd": set(), "rev": set()}
    combos = [(i, j) for i in range(96) for j in range(96)
              if i != j and i not in exclude["fwd"] and j not in exclude["rev"]]
    if len(combos) < len(wells):
        raise ValueError(f"{len(wells)} wells but only {len(combos)} barcode combinations")
    order = rng.permutation(len(combos))[:len(wells)]
    return [(p, w, combos[k][0], combos[k][1]) for (p, w), k in zip(wells, order)]


def minibar_sheet(assigned, table):
    fwd = table[table.side == "fwd"].set_index("idx")
    rev = table[table.side == "rev"].set_index("idx")
    return pd.DataFrame([dict(SampleID=f"Plate{int(p):02d}_{w}",
                              FwIndex=fwd.at[i, "barcode"], FwPrimer=FWD_PRIMER,
                              RvIndex=rev.at[j, "barcode"], RwPrimer=REV_PRIMER)
                         for p, w, i, j in assigned])


def min_edit_distance(seqs):
    import edlib
    seqs = sorted(set(seqs))
    return min(edlib.align(a, b, task="distance")["editDistance"]
               for k, a in enumerate(seqs) for b in seqs[k + 1:])


def demux_config(barcodes_tsv, outdir, fastq_input="<fill in after the run>",
                 barcode_edit_distance=4):
    """Written from the same inputs as the sheet, with the key names demultiplex.smk
    actually reads (`barcode_edit_distance`, not `barcode_edit_dist`, which was silently
    ignored in every earlier config)."""
    return {
        "fastq_input": fastq_input,
        "barcodes_tsv": str(barcodes_tsv),
        "nanofilt": {"min_len": 1400, "max_len": 1800, "min_q": 10, "headcrop": 0, "tailcrop": 0},
        "minibar": {"prefix": "", "barcode_edit_distance": barcode_edit_distance,
                    "primer_edit_distance": 6, "search_length": 150},
        "outdir": str(outdir),
    }
