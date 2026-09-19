process SANITIZEME {

    tag "${reads.simpleName}"

    cpus 8

    publishDir "${params.outdir}/03_removehost",
        mode: 'copy'

    input:

    path reads

    path human_reference

    output:

    path "${reads.simpleName}_filtered.fastq"

    script:

    """
    SanitizeMe_CLI.py \
        -i . \
        -r ${human_reference} \
        -o . \
        -t ${task.cpus} \
        --Nanopore
    """
}
