process DORADO_DEMUX {

    tag "demultiplex"

    publishDir "${params.outdir}/02_demultiplexed", mode: 'copy'

    input:
    path bam

    output:
    path "demultiplexed/**/*.fastq"

    script:
    """
    dorado demux \
        ${bam} \
        --emit-fastq \
        --kit-name ${params.kit_name} \
        --output-dir demultiplexed
    """
}
