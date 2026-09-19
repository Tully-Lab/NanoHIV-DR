process HIVDB_READS {

    tag "${sample}"

    cpus 1

    publishDir "${params.outdir}/08_hivdb/reads",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample), path(codfreq)

    output:
    tuple val(sample), path("${sample}.reads.json"),
        emit: reports

    script:
    """
    set -euo pipefail

    echo "============================================================"
    echo " STANFORD HIVDB BY-READS"
    echo "============================================================"
    echo ""
    echo "Sample: ${sample}"
    echo "CodFreq: ${codfreq}"
    echo ""
    echo "Parameters:"
    echo "  Minimum prevalence:       10%"
    echo "  Maximum mixture rate:      2%"
    echo "  Minimum codon reads:       1"
    echo "  Minimum position reads:   50"
    echo ""

    cp -L "${codfreq}" input.codfreq

    micromamba run -n report sierrapy seqreads \
        -p 0.10 \
        -m 0.02 \
        -d 1 \
        -D 50 \
        input.codfreq

    if [[ ! -s input.report.json ]]; then
        echo "ERROR: SierraPy did not produce input.report.json"
        ls -lah
        exit 1
    fi

    micromamba run -n report python -m json.tool \
        input.report.json >/dev/null

    mv input.report.json "${sample}.reads.json"

    echo ""
    echo "HIVdb By-Reads result:"
    ls -lh "${sample}.reads.json"
    echo ""
    """
}
