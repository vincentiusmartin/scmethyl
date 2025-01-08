
process METHYLDACKEL {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/methyldackel", mode: 'copy'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Vincentius/envs/biscuit'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.bedGraph"), emit: meth
    tuple val(meta), path("*.svg"), optional: true, emit: mbias


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -f  ${prefix}.bam ] && ln -sf ${read} ${prefix}.bam
    
    samtools index ${prefix}.bam
    MethylDackel mbias ${params.genome_ref} ${prefix}.bam ${prefix}
    MethylDackel extract -o "${prefix}" --mergeContext --minDepth 1 --OT 10,145,10,145 --OB 10,145,10,145 ${params.genome_ref} ${prefix}.bam
    """
    
}