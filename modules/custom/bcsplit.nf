
process BCSPLIT {
  tag "$meta.id"
  label 'process_high'
  
  module 'mamba'
  conda '/research/groups/northcgrp/home/common/Vincentius/envs/renv'
  
  publishDir "${params.outdir}/custom/bcsplit", mode: 'copy'
  
  input:
  tuple val(meta), path(reads)
    
  output:
  tuple val(meta), path("*.fastq.gz"), emit: reads
  tuple val(meta), path("*.bcsplit.txt"), emit: log
    
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
  
  def read1 = []
  readList.eachWithIndex { v, ix -> if (!(ix & 1)) read1 << v }
  
  """
  bcsplit.R ${read1.join(',')} $prefix
  """
}