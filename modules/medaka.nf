/*
 * ============================================================
 * MEDAKA
 * ============================================================
 *
 * Input:
 *
 *     NanoQ filtered FASTQ
 *             +
 *     Reference sequence
 *
 * Output:
 *
 *     consensus.fasta
 *
 * ============================================================
 */

process MEDAKA {

    tag "${reads.simpleName}"

    cpus params.threads

    publishDir "${params.outdir}/05_medaka",
        mode: 'copy',
        overwrite: true

    input:

    tuple path(reference), path(reads)


    output:

    tuple val(reads.simpleName), path("consensus.fasta"),
        emit: consensus


    script:

    """
    set -euo pipefail

    echo "============================================================"
    echo " MEDAKA"
    echo "============================================================"
    echo ""

    echo "Sample:"
    echo "    ${reads.simpleName}"
    echo ""

    echo "Input FASTQ:"
    echo "    ${reads}"
    echo ""

    echo "Reference:"
    echo "    ${reference}"
    echo ""

    echo "Medaka model:"
    echo "    ${params.medaka_model}"
    echo ""

    echo "Threads:"
    echo "    ${task.cpus}"
    echo ""

    # --------------------------------------------------------
    # Run Medaka
    # --------------------------------------------------------
    #
    # medaka_consensus accepts the filtered FASTQ directly.
    #
    # --------------------------------------------------------

    medaka_consensus \
        -i "${reads}" \
        -d "${reference}" \
        -o medaka \
        -t ${task.cpus} \
        -m "${params.medaka_model}"


    # --------------------------------------------------------
    # Validate Medaka output
    # --------------------------------------------------------

    if [[ ! -f medaka/consensus.fasta ]]; then

        echo ""
        echo "ERROR: Medaka did not produce:"
        echo ""
        echo "    medaka/consensus.fasta"
        echo ""

        echo "Contents of Medaka output directory:"
        ls -lah medaka 2>/dev/null || true

        exit 1

    fi


    if [[ ! -s medaka/consensus.fasta ]]; then

        echo ""
        echo "ERROR: Medaka consensus.fasta is empty."
        echo ""

        exit 1

    fi


    # --------------------------------------------------------
    # Create standard pipeline output
    # --------------------------------------------------------

    cp \
        medaka/consensus.fasta \
        consensus.fasta


    # --------------------------------------------------------
    # Validate FASTA
    # --------------------------------------------------------

    if ! grep -q '^>' consensus.fasta; then

        echo ""
        echo "ERROR: Medaka consensus is not a valid FASTA file."
        echo ""

        exit 1

    fi


    # --------------------------------------------------------
    # Report output
    # --------------------------------------------------------

    echo ""
    echo "============================================================"
    echo " MEDAKA COMPLETE"
    echo "============================================================"
    echo ""

    echo "Consensus:"
    echo "    consensus.fasta"
    echo ""

    ls -lh consensus.fasta

    echo ""
    """
}
