# NanoHIV-DR reporting

This folder holds the clinical reporting component of NanoHIV-DR. It is invoked
automatically by the Nextflow `REPORT` process in Stage 2 (`--stage report`) and
runs inside the pipeline's Docker container, so it does not normally need to be
run by hand.

It is adapted from Samantha Campbell's UVRI-HIV-diagnostic-report workflow
(https://github.com/centre-for-virus-research/UVRI-HIV-diagnostic-report). The
original phylogeny step (MAFFT alignment + RAxML tree) has been **removed** —
NanoHIV-DR does not perform phylogenetic analysis.

## What it does

`preprocessing.sh` takes the Stage 2 inputs and generates the per-sample
clinical reports:

1. Collects the Stanford HIVdb results — the consensus analysis JSON and the
   per-sample by-reads JSON files.
2. Runs `bin/parse_json_write_docx.py` to produce a Word (`.docx`) clinical
   report per sample and a subtype table.

## Usage

```
preprocessing.sh \
    -f <consensus.fasta> \
    -t <metadata.tsv> \
    -j <hivdb_consensus.json> \
    -r <hivdb_reads_directory>
```

| Flag | Description |
| --- | --- |
| `-f` | Consensus FASTA |
| `-t` | NanoHIV-DR metadata TSV |
| `-j` | Stanford HIVdb consensus JSON |
| `-r` | Directory of per-sample HIVdb by-reads JSON (`*.reads.json`) |
| `-h` | Show help |

## Output

All results are written to `results_report/`:

- `clinical_reports/` — per-sample clinical reports (`.docx`)
- `subtypes.tsv` — HIV-1 subtype per sample
- `hivdb_results.json` — Stanford HIVdb consensus analysis
- `hivdb_reads/` — Stanford HIVdb by-reads analyses
- `consensus.fasta`, `metadata.tsv` — copies of the inputs used

## Requirements

Provided by the pipeline's Docker image (see the top-level `Dockerfile`):

- Python 3.10 with `sierra-client` (SierraPy), `python-docx`, and `pandas`
  (see `requirements.txt`)
