process REPORT {

    tag "clinical-reports"

    cpus 4

    publishDir "${params.outdir}/09_reports",
        mode: 'copy',
        overwrite: true

    input:

    path fasta
    path metadata
    path consensus_json
    path reads_json

    path codfreq_files

    output:

    path "results_report"

    script:

    """
    set -euo pipefail

    mkdir -p hivdb_reads
    cp ${reads_json} hivdb_reads/

    mkdir -p codfreq_reads
    cp ${codfreq_files} codfreq_reads/ 2>/dev/null || true

    echo "============================================================"
    echo " REPORT INPUTS"
    echo "============================================================"
    echo "Consensus:       ${fasta}"
    echo "Metadata:        ${metadata}"
    echo "Consensus HIVdb: ${consensus_json}"
    echo "By-Reads HIVdb:"
    ls -1 hivdb_reads/*.reads.json
    echo ""

    /opt/REPORT/preprocessing.sh \
        -f "${fasta}" \
        -t "${metadata}" \
        -j "${consensus_json}" \
        -r hivdb_reads \
        -c codfreq_reads
    """
}