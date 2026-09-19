/*
 * ============================================================
 * MINIMAP2
 * ============================================================
 *
 * Input:
 *
 *     NanoQ filtered FASTQ
 *             +
 *     HIV-1 minimap2 index (.mmi)
 *
 * Output:
 *
 *     BAM
 *     BAI
 *
 * This is Branch 2 of the NanoHIV-DR workflow.
 *
 * NanoQ
 *    ↓
 * MINIMAP2
 *    ↓
 * BAM + BAI
 *
 * ============================================================
 */

process MINIMAP2 {

    tag "${reads.simpleName}"

    cpus params.threads

    publishDir "${params.outdir}/06_minimap2",
        mode: 'copy',
        overwrite: true

    input:

    tuple path(index), path(reads)

    output:

    tuple path("*.bam"), path("*.bam.bai"),
        emit: alignment

    script:

    def sample = reads.simpleName
        .replaceFirst(/\.trimmed$/, '')

    """
    set -euo pipefail

    echo "============================================================"
    echo " MINIMAP2"
    echo "============================================================"
    echo ""

    echo "Sample:"
    echo "    ${sample}"
    echo ""

    echo "Input FASTQ:"
    echo "    ${reads}"
    echo ""

    echo "Reference index:"
    echo "    ${index}"
    echo ""

    echo "Threads:"
    echo "    ${task.cpus}"
    echo ""

    # --------------------------------------------------------
    # Alignment
    # --------------------------------------------------------

    minimap2 \
        -ax map-ont \
        -t ${task.cpus} \
        "${index}" \
        "${reads}" \
        | samtools sort \
            -@ ${task.cpus} \
            -o "${sample}.bam" \
            -

    # --------------------------------------------------------
    # BAM index
    # --------------------------------------------------------

    samtools index \
        -@ ${task.cpus} \
        "${sample}.bam"

    # --------------------------------------------------------
    # Validate outputs
    # --------------------------------------------------------

    if [[ ! -s "${sample}.bam" ]]; then

        echo ""
        echo "ERROR: Minimap2 produced an empty BAM file."
        echo ""

        exit 1

    fi

    if [[ ! -s "${sample}.bam.bai" ]]; then

        echo ""
        echo "ERROR: Samtools failed to create BAM index."
        echo ""

        exit 1

    fi

    # --------------------------------------------------------
    # BAM statistics
    # --------------------------------------------------------

    echo ""
    echo "BAM statistics:"
    echo ""

    samtools flagstat \
        "${sample}.bam" \
        || true

    # --------------------------------------------------------
    # Final output
    # --------------------------------------------------------

    echo ""
    echo "============================================================"
    echo " MINIMAP2 COMPLETE"
    echo "============================================================"
    echo ""

    echo "BAM:"
    echo "    ${sample}.bam"
    echo ""

    echo "BAI:"
    echo "    ${sample}.bam.bai"
    echo ""

    ls -lh \
        "${sample}.bam" \
        "${sample}.bam.bai"

    echo ""
    """
}
