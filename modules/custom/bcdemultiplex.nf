

process BCDEMUL {
  tag "$meta.id"
  label 'process_high'
  
  module 'mamba'
  conda '/research/groups/northcgrp/home/common/Vincentius/envs/scmethyl'
  
  publishDir "${params.outdir}/custom/bcdemul", mode: 'copy'
  
  input:
  tuple val(meta), path(reads)
    
  output:
  tuple val(meta), path("*.fastq.gz"), emit: reads
  tuple val(meta), path("*_counts.csv"), emit: log
    
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
  
  def read1 = []
  def read2 = []
  readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
  
  """
  bcdemultiplex.py --fq1 ${read1.join(' ')} --fq2 ${read2.join(' ')} \
    --sample_name $prefix --bc ${params.barcode_path}
  """
}