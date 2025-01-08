
process SAMBAMBA_MARKDUP {
  tag "$meta.id"
  label 'process_med'
  
  module 'mamba'
  conda '/research/groups/northcgrp/home/common/Vincentius/envs/methylctools'
  
  publishDir "${params.outdir}/sambamba/markdup", mode: 'copy'
  
  input:
  tuple val(meta), path(read)
  
  output:
  tuple val(meta), path("*.markdup.bam"), emit: reads
  tuple val(meta), path("*.markdup.bam.flagstat"), emit: flagstat

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"

  """
  [ ! -f  ${prefix}.bam ] && ln -sf ${read} ${prefix}.bam
  
  samtools index ${prefix}.bam
  sambamba markdup ${prefix}.bam ${prefix}.markdup.bam
  samtools flagstat ${prefix}.markdup.bam > ${prefix}.markdup.bam.flagstat
  """
}