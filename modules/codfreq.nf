/*
 * ============================================================
 * CODFREQ
 * ============================================================
 *
 * Input:
 *
 *     NanoQ filtered FASTQ
 *             +
 *     CodFreq profile (.json)
 *
 * Output:
 *
 *     *.codfreq
 *
 * ============================================================
 */

process CODFREQ {

    tag "${reads.simpleName}"

    cpus params.threads

    publishDir "${params.outdir}/07_codfreq",
        mode: 'copy',
        overwrite: true

    input:
    tuple path(reads), path(profile)

    output:
    tuple val(reads.simpleName), path("*.codfreq"),
        emit: codfreq

    script:

    def sample = reads.simpleName
        .replaceFirst(/\.trimmed$/, '')

    """
    set -euo pipefail

    echo "============================================================"
    echo " CODFREQ"
    echo "============================================================"
    echo ""

    echo "Sample:"
    echo "    ${sample}"
    echo ""

    echo "Input FASTQ:"
    echo "    ${reads}"
    echo ""

    echo "CodFreq profile:"
    echo "    ${profile}"
    echo ""

    echo "Threads:"
    echo "    ${task.cpus}"
    echo ""

    # --------------------------------------------------------
    # Check CodFreq installation
    # --------------------------------------------------------

    if ! command -v fastq2codfreq >/dev/null 2>&1; then
        echo ""
        echo "ERROR: fastq2codfreq was not found in PATH."
        echo ""
        echo "PATH:"
        echo "\${PATH}"
        exit 1
    fi

    echo "CodFreq executable:"
    command -v fastq2codfreq
    echo ""

    # --------------------------------------------------------
    # Create CodFreq working directory
    # --------------------------------------------------------

    mkdir -p codfreq_work

    cp "${reads}" codfreq_work/

    # --------------------------------------------------------
    # Run CodFreq
    # --------------------------------------------------------

    fastq2codfreq \
        --program minimap2 \
        --profile "${profile}" \
        --workers ${task.cpus} \
        codfreq_work

    # --------------------------------------------------------
    # Locate and validate output
    # --------------------------------------------------------

    echo ""
    echo "CodFreq output files:"
    echo ""

    find codfreq_work -maxdepth 1 -type f -print

    codfreq_file=\$(find codfreq_work \
        -maxdepth 1 \
        -type f \
        -name "*.codfreq" \
        | head -n 1)

    if [[ -z "\${codfreq_file}" ]]; then
        echo ""
        echo "ERROR: No *.codfreq output was produced."
        echo ""
        find codfreq_work -maxdepth 1 -type f -print
        exit 1
    fi

    # --------------------------------------------------------
    # Create standard pipeline output
    # --------------------------------------------------------

    cp "\${codfreq_file}" "${sample}.codfreq"

    if [[ ! -s "${sample}.codfreq" ]]; then
        echo ""
        echo "ERROR: CodFreq result is empty."
        echo ""
        exit 1
    fi

    echo ""
    echo "============================================================"
    echo " CODFREQ COMPLETE"
    echo "============================================================"
    echo ""

    echo "CodFreq result:"
    echo "    ${sample}.codfreq"
    echo ""

    ls -lh "${sample}.codfreq"
    echo ""
    """
}