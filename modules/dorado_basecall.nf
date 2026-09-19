process DORADO_BASECALL {

    tag "basecalling"

    cpus 8

    publishDir "${params.outdir}/01_basecalled", mode: 'copy'

    input:
    path pod5_dir
    path dorado_model

    output:
    path "all_reads.bam"

    script:
    """
    dorado basecaller \
        ${dorado_model} \
        ${pod5_dir} \
        --device ${params.device} \
        --kit-name ${params.kit_name} \
        > all_reads.bam
    """
}