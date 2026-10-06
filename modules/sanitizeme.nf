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
    # SanitizeMe is installed in its own micromamba environment ('host')
    # so that its Gooey/wxPython GUI dependency stack does not have to
    # co-solve with the base bioinformatics tools (minimap2, samtools,
    # ivar, ...). Installed together they make the aarch64 image build
    # unsolvable; isolating SanitizeMe keeps it available unchanged.
    micromamba run -n host SanitizeMe_CLI.py \
        -i . \
        -r ${human_reference} \
        -o . \
        -t ${task.cpus} \
        --Nanopore
    """
}
