
process METHYLCTOOLS_FQCONV {
  tag "$meta.id"
  label 'process_med'
  
  module 'samtools/1.10'
  module 'bwa/0.7.17'
  module 'mamba'
  conda '/research/groups/northcgrp/home/common/Vincentius/envs/scWGBS_pip'
  
  publishDir "${params.outdir}/methylctools/fqconv", mode: 'copy'
  
  input:
  tuple val(meta), path(reads)
  
  output:
  tuple val(meta), path("*.bam")
  //tuple val(meta), path("*metrics.txt"), emit: metrics
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"

  """
  [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
  [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
  
  methylCtools fqconv -1 ${prefix}_2.fastq.gz -2 ${prefix}_1.fastq.gz - | \
    bwa mem -t 4 -T 0 -M -p ${params.genome_ref_conv} - | \
    samtools view -Sbu - | \
    methylCtools bconv -m ${prefix}.metrics.txt - - | \
    samtools sort - > ${prefix}.bam
  """
  // 
}