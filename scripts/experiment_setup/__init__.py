"""Shared, importable experiment-setup library (see ../../experiments/README.md).

Every per-experiment input lives in that experiment's `01_setup/setup_config.yaml`; its
`run_setup.py` calls into here, writes `outputs/`, and runs `validate`. Nothing in this
package writes a file on import or holds experiment-specific constants.

    plates       plate geometry: usable wells, neighbours
    pool         which strains may be used, one well per distinct genome
    communities  community generators (size x relatedness, series, tests)
    layout       communities x replicates -> plates and wells
    barcodes     wells -> forward x reverse barcode combinations -> minibar sheet
    echo         transfer plan -> Echo CSVs, with OD-normalised volumes and refill batches
    validate     checks on every output; run_setup refuses to finish if one fails
"""
