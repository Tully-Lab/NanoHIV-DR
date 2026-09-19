/*
 * ============================================================
 * IVAR CONSENSUS
 * ============================================================
 *
 * Minimap2 BAM + BAI
 *        +
 * HXB2 reference
 *        ↓
 * samtools mpileup
 *        ↓
 * iVar consensus
 *
 * ============================================================
 */

process IVAR_CONSENSUS {

    tag "ivar-consensus"

    publishDir "${params.outdir}/05_consensus",
        mode: 'copy',
        overwrite: true

    input:

    tuple path(bam), path(bai), path(reference)

    output:

    path "*.consensus.fa",
        emit: consensus

    script:

    sample = bam.simpleName
        .replaceFirst(/_filtered$/, '')

    """
    set -euo pipefail

    echo "============================================================"
    echo " IVAR CONSENSUS"
    echo "============================================================"
    echo ""

    echo "Sample:"
    echo "    ${sample}"
    echo ""

    echo "Input BAM:"
    echo "    ${bam}"
    echo ""

    echo "Reference:"
    echo "    ${reference}"
    echo ""

    echo "Parameters:"
    echo "    Minimum base quality: ${params.ivar_min_quality}"
    echo "    Consensus threshold:  ${params.ivar_threshold}"
    echo "    Minimum depth:        ${params.ivar_min_depth}"
    echo ""

    samtools mpileup \
        -aa \
        -A \
        -d ${params.ivar_max_depth} \
        -Q 0 \
        -f "${reference}" \
        "${bam}" \
        | ivar consensus \
            -p "${sample}" \
            -q ${params.ivar_min_quality} \
            -t ${params.ivar_threshold} \
            -m ${params.ivar_min_depth} \
            -n N

    if [[ ! -s "${sample}.fa" ]]; then
        echo ""
        echo "ERROR: iVar did not produce ${sample}.fa"
        echo ""
        ls -lah
        exit 1
    fi

	# --------------------------------------------------------
	# Normalise FASTA header
	# --------------------------------------------------------

	sed "1s/.*/>${sample}/" \
    	"${sample}.fa" \
    	> "${sample}.consensus.fa"

	rm "${sample}.fa"
    echo ""
    echo "Consensus statistics:"
    echo ""

    awk '
        !/^>/ {
            seq = seq \$0
        }
        END {
            tmp = seq
            n = gsub(/[Nn]/, "", tmp)
            callable = length(seq) - n

            print "    Length:   " length(seq)
            print "    Callable: " callable
            print "    Ns:       " n

            if (length(seq) > 0) {
                printf "    Callable: %.1f%%\\n", \
                    100 * callable / length(seq)
            }
        }
    ' "${sample}.consensus.fa"

    echo ""
    echo "============================================================"
    echo " IVAR CONSENSUS COMPLETE"
    echo "============================================================"
    echo ""

    ls -lh "${sample}.consensus.fa"
    echo ""
    """
}
