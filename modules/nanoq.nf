process NANOQ {

    tag "${reads.simpleName}"

    cpus 8

    publishDir "${params.outdir}/04_nanoq",
        mode: 'copy',
        overwrite: true

    input:

    path reads

    output:

    path "${reads.simpleName}.trimmed.fastq"

    script:

    """
    set -euo pipefail

    echo "============================================================"
    echo " NanoQ"
    echo "============================================================"

    echo "Input:       ${reads}"
    echo "Min quality: ${params.min_quality}"
    echo "Min length:  ${params.min_length}"
    echo "Max length:  ${params.max_length}"
    echo ""

    nanoq \\
        --input "${reads}" \\
        --min-qual ${params.min_quality} \\
        --min-len ${params.min_length} \\
        --max-len ${params.max_length} \\
        --output "${reads.simpleName}.trimmed.fastq"

    if [[ ! -s "${reads.simpleName}.trimmed.fastq" ]]; then

        echo ""
        echo "ERROR: NanoQ produced an empty output file."
        echo ""

        exit 1

    fi

    echo ""
    echo "NanoQ filtering complete."
    echo "Output: ${reads.simpleName}.trimmed.fastq"
    echo ""
    """
}
