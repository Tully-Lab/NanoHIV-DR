/*
 * ============================================================
 * DRM-LINK  (read-backed DRM co-occurrence / within-amplicon linkage)
 * ============================================================
 *
 * EXPERIMENTAL / OPTIONAL. Not wired into the default workflow.
 * Consumes the per-sample coordinate-sorted BAM (from MINIMAP2) and a
 * DRM target list, and reports pairwise within-amplicon linkage.
 *
 * The tool (bin/drm_link.py) is placed on the task PATH automatically by
 * Nextflow. It needs pysam in the container's python (added to the base
 * environment in the Dockerfile). See bin/drm_link.py for the method and
 * the DRM-LINK spec for the claim language.
 *
 * To enable, in main.nf Stage 1 after MINIMAP2:
 *
 *   targets_ch = channel.value(file(params.drm_targets, checkIfExists: true))
 *   LINKAGE(alignment.alignment.map { bam, bai ->
 *               tuple(bam.simpleName.replaceFirst(/_filtered$/,''), bam, bai) },
 *           targets_ch)
 * ============================================================
 */

process LINKAGE {

    tag "${sample}"

    publishDir "${params.outdir}/10_linkage",
        mode: 'copy',
        overwrite: true

    input:

    tuple val(sample), path(bam), path(bai)

    path targets

    output:

    path "linkage_out/*.linkage.tsv",
        emit: linkage

    script:

    """
    drm_link.py \
        --bam ${bam} \
        --targets ${targets} \
        --sample ${sample} \
        --outdir linkage_out \
        --chimera-floor ${params.chimera_floor} \
        --min-mapq ${params.linkage_min_mapq} \
        --min-cocov ${params.linkage_min_cocov}
    """
}
