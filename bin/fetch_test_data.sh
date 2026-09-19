#!/usr/bin/env bash
#
# Fetch a small public HIV-1 Nanopore test set into data/test/.
#
# Source: Ito et al. 2024, near-full-length HIV-1 Nanopore sequencing
#         (DDBJ DRA BioProject PRJDB17699). These reads span the whole
#         genome, so they cover the full pol gene (PR + RT + IN).
#
# Two download methods are tried, in order:
#   1. SRA toolkit  (prefetch + fasterq-dump)   <- most reliable for DRR
#   2. ENA FASTQ mirror (curl)                   <- fallback
#
# Install the SRA toolkit if you don't have it:
#     conda install -c bioconda sra-tools
#
# Override the accession list as arguments, e.g.:
#     bin/fetch_test_data.sh DRR537715 DRR537716 DRR537717
#
set -euo pipefail

OUTDIR="data/test"
mkdir -p "$OUTDIR"

ACCESSIONS=("$@")
if [[ ${#ACCESSIONS[@]} -eq 0 ]]; then
    ACCESSIONS=(DRR537715 DRR537716)
fi

fetch_sra() {
    local acc="$1"
    echo "  [SRA toolkit] prefetch + fasterq-dump ${acc} ..."
    prefetch "$acc" -O "$OUTDIR" >/dev/null
    fasterq-dump "$acc" -O "$OUTDIR" >/dev/null
    # fasterq-dump writes ${acc}.fastq (single-end / Nanopore) -> gzip it
    for f in "$OUTDIR/${acc}".fastq "$OUTDIR/${acc}"_*.fastq; do
        [[ -f "$f" ]] && gzip -f "$f"
    done
    rm -rf "$OUTDIR/${acc}"    # remove the prefetch .sra cache dir
}

fetch_ena() {
    local acc="$1"
    # Ask for both columns; the FASTQ URLs live in the fastq_ftp column ($2)
    local ftp
    ftp=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp" \
          | awk -F'\t' 'NR==2{print $2}')
    if [[ -z "${ftp:-}" ]]; then
        return 1
    fi
    local IFS=';'
    for u in $ftp; do
        [[ -z "$u" ]] && continue
        echo "  [ENA] downloading https://${u}"
        wget -q --show-progress -P "$OUTDIR" "https://${u}"
    done
}

for acc in "${ACCESSIONS[@]}"; do
    echo ""
    echo "== ${acc} =="
    if command -v fasterq-dump >/dev/null 2>&1 && command -v prefetch >/dev/null 2>&1; then
        fetch_sra "$acc" || { echo "  SRA toolkit failed; trying ENA ..."; fetch_ena "$acc" || echo "  ERROR: could not fetch ${acc} from SRA or ENA."; }
    else
        echo "  SRA toolkit not found; trying ENA mirror ..."
        fetch_ena "$acc" || {
            echo "  ERROR: could not fetch ${acc} from ENA, and the SRA toolkit is not installed."
            echo "  Install it and re-run:  conda install -c bioconda sra-tools"
        }
    fi
done

echo ""
echo "Files in ${OUTDIR}/:"
ls -lh "$OUTDIR"/*.fastq* 2>/dev/null || echo "  (none downloaded — see messages above)"
echo ""
echo "If files are present, run:  nextflow run main.nf -profile test"
