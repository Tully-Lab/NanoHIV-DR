#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REPORT_PYTHON=(micromamba run -n report python)

FASTA=""
INFO=""
HIVDB_JSON=""
READS_DIR=""
CODFREQ_DIR=""

usage() {

    cat <<USAGE
Usage:

    preprocessing.sh \
        -f <consensus.fasta> \
        -t <metadata.tsv> \
        -j <hivdb_consensus.json> \
        -r <hivdb_reads_directory>

Required:

    -f    Consensus FASTA
    -t    NanoHIV-DR metadata TSV
    -j    Stanford HIVdb consensus JSON
    -r    Directory containing per-sample HIVdb By-Reads JSON files
    -h    Show this help

USAGE
}

while getopts ":hf:t:j:r:c:" opt; do

    case "$opt" in

        h)
            usage
            exit 0
            ;;

        f)
            FASTA="$OPTARG"
            ;;

        t)
            INFO="$OPTARG"
            ;;

        j)
            HIVDB_JSON="$OPTARG"
            ;;

        r)
            READS_DIR="$OPTARG"
            ;;

        c)
            CODFREQ_DIR="$OPTARG"
            ;;

        :)
            echo "ERROR: Option -$OPTARG requires an argument." >&2
            exit 1
            ;;

        \?)
            echo "ERROR: Invalid option -$OPTARG" >&2
            usage
            exit 1
            ;;

    esac

done


# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

if [[ -z "$FASTA" ]]; then
    echo "ERROR: Missing FASTA (-f)." >&2
    exit 1
fi

if [[ -z "$INFO" ]]; then
    echo "ERROR: Missing metadata TSV (-t)." >&2
    exit 1
fi

if [[ -z "$HIVDB_JSON" ]]; then
    echo "ERROR: Missing consensus HIVdb JSON (-j)." >&2
    exit 1
fi

if [[ -z "$READS_DIR" ]]; then
    echo "ERROR: Missing HIVdb By-Reads directory (-r)." >&2
    exit 1
fi


if [[ ! -s "$FASTA" ]]; then
    echo "ERROR: FASTA not found or empty: $FASTA" >&2
    exit 1
fi

if [[ ! -s "$INFO" ]]; then
    echo "ERROR: Metadata not found or empty: $INFO" >&2
    exit 1
fi

if [[ ! -s "$HIVDB_JSON" ]]; then
    echo "ERROR: Consensus HIVdb JSON not found or empty: $HIVDB_JSON" >&2
    exit 1
fi

if [[ ! -d "$READS_DIR" ]]; then
    echo "ERROR: HIVdb By-Reads directory not found: $READS_DIR" >&2
    exit 1
fi


NSEQ=$(grep -c '^>' "$FASTA" || true)

if [[ "$NSEQ" -lt 1 ]]; then
    echo "ERROR: No FASTA sequences found in $FASTA" >&2
    exit 1
fi


NREADS=$(find "$READS_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.reads.json' \
    | wc -l | tr -d ' ')

if [[ "$NREADS" -lt 1 ]]; then
    echo "ERROR: No *.reads.json files found in $READS_DIR" >&2
    exit 1
fi


# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

RESULTS="results_report"
REPORT_DIR="clinical_reports"

JSON="${RESULTS}/hivdb_results.json"
READS_RESULTS="${RESULTS}/hivdb_reads"
SUBTYPES="${RESULTS}/subtypes.tsv"

rm -rf "$RESULTS" "$REPORT_DIR"

mkdir -p "$RESULTS"
mkdir -p "$READS_RESULTS"


echo "============================================================"
echo " NanoHIV-DR clinical reporting"
echo "============================================================"
echo ""
echo "Input FASTA:          $FASTA"
echo "Metadata:             $INFO"
echo "Sequences:            $NSEQ"
echo "Consensus HIVdb JSON: $HIVDB_JSON"
echo "By-Reads results:     $NREADS"
echo ""


# ------------------------------------------------------------
# Preserve Stanford HIVdb outputs
# ------------------------------------------------------------

echo "Collecting Stanford HIVdb results..."

cp "$HIVDB_JSON" "$JSON"
cp "$READS_DIR"/*.reads.json "$READS_RESULTS"/

echo ""
echo "Consensus HIVdb result:"
echo "    $JSON"

echo ""
echo "By-Reads HIVdb results:"
ls -1 "$READS_RESULTS"/*.reads.json

echo ""


# ------------------------------------------------------------
# Generate clinical reports
# ------------------------------------------------------------

echo "Generating clinical reports..."

"${REPORT_PYTHON[@]}" \
    "${script_dir}/bin/parse_json_write_docx.py" \
    --json "$JSON" \
    --reads-dir "$READS_RESULTS" \
    --data "$INFO" \
    --fasta "$FASTA" \
    --output "$SUBTYPES" \
    --reports \
    --report-dir "$REPORT_DIR" \
    ${CODFREQ_DIR:+--codfreq-dir "$CODFREQ_DIR"}

if [[ -d "$REPORT_DIR" ]]; then
    mv "$REPORT_DIR" "$RESULTS/"
fi


# ------------------------------------------------------------
# Preserve source inputs
# ------------------------------------------------------------

cp "$FASTA" "${RESULTS}/consensus.fasta"
cp "$INFO" "${RESULTS}/metadata.tsv"


echo ""
echo "============================================================"
echo " NanoHIV-DR reporting complete"
echo "============================================================"
echo ""
echo "Consensus HIVdb JSON:"
echo "    $JSON"
echo ""
echo "By-Reads HIVdb JSON:"
echo "    $READS_RESULTS/"
echo ""
echo "Subtype table:"
echo "    $SUBTYPES"
echo ""
echo "Clinical reports:"
echo "    ${RESULTS}/${REPORT_DIR}/"
echo ""
echo "============================================================"