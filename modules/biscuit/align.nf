/**
 * Perform biscuit aligner and convertion to sorted bam file
 */
process BISCUIT_ALIGN {
    tag "$meta.id"
    label 'process_med'

    module 'mamba'
    conda '/research/groups/northcgrp/home/common/Vincentius/envs/biscuit'

    publishDir "${params.outdir}/biscuit/align", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    
    biscuit align -@ ${task.cpus} -R "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA\\tLB:LIB" \
        ${params.genome_ref} \
        ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz | \
    samtools view -@ ${task.cpus} -S -b |
    samtools sort -@ ${task.cpus} -o ${prefix}.bam 
    """

}
