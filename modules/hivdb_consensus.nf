process HIVDB_CONSENSUS {

    tag "consensus"

    cpus 1

    publishDir "${params.outdir}/08_hivdb/consensus",
        mode: 'copy',
        overwrite: true

    input:
    path fasta

    output:
    path "hivdb_consensus.json",
        emit: json

    script:
    """
    set -euo pipefail

    echo "============================================================"
    echo " STANFORD HIVDB CONSENSUS"
    echo "============================================================"
    echo ""
    echo "Input FASTA: ${fasta}"
    echo ""

    micromamba run -n report sierrapy fasta \
        "${fasta}" \
        -o hivdb_consensus.json \
        --no-sharding

    if [[ ! -s hivdb_consensus.json ]]; then
        echo "ERROR: SierraPy did not produce consensus HIVdb JSON."
        exit 1
    fi

    micromamba run -n report python -m json.tool \
        hivdb_consensus.json >/dev/null

    echo ""
    echo "Consensus HIVdb result:"
    ls -lh hivdb_consensus.json
    echo ""
    """
}
