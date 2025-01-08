
process METHYLCTOOLS_BCALL {
  tag "$meta.id"
  label 'process_med'
  
  module 'mamba'
  conda '/research/groups/northcgrp/home/common/Vincentius/envs/methylctools'
  
  publishDir "${params.outdir}/methylctools/bcall", mode: 'copy'
  
  input:
  tuple val(meta), path(read)
  
  output:
  tuple val(meta), path("*")

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"

  """
  [ ! -f  ${prefix}.bam ] && ln -sf ${read} ${prefix}.bam
  
  methylCtools bcall --silent --trimPE -m ${prefix}.call.metrics ${params.genome_cpos} ${prefix}.bam - | \
    gzip > ${prefix}.call.gz
  methylCtools bcall --silent --trimPE -r chr1:1-10000000 -m ${prefix}.callCH.metrics ${params.genome_cpos} ${prefix}.bam - | \
    gzip > ${prefix}.callCH.gz
  """
}