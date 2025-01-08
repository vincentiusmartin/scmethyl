
process FASTP {
  tag "$meta.id"
  label 'process_med'
  
  module 'fastp/0.19.6'
  
  publishDir "${params.outdir}/fastp/fastp", mode: 'copy'
  
  input:
  tuple val(meta), path(reads)
    
  output:
  tuple val(meta), path("*.fastp.fastq.gz"), emit: reads
  tuple val(meta), path("*fastp.html"), emit: html
  tuple val(meta), path("*fastp.json"), emit: json
    
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  
  def out_fq1 = meta.single_end ? "${prefix}.fastp.fastq.gz" : "${prefix}_1.fastp.fastq.gz" 
  def out_fq2 = "${prefix}_2.fastp.fastq.gz"
  
  def args_list = ["--adapter_sequence=${params.adapter_sequence}",
    "--adapter_sequence_r2=${params.adapter_sequence_r2}", 
    "--trim_front1 ${params.trim_front}", "--trim_front2 ${params.trim_front}",
    "--trim_tail1 ${params.trim_tail}", "--trim_tail2 ${params.trim_tail}",
    "--trim_poly_g"]
  
  if (!meta.single_end) {
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    fastp ${args_list.join(' ')} \
      -i ${prefix}_1.fastq.gz -I ${prefix}_2.fastq.gz \
      -o $out_fq1 -O $out_fq2 \
      --html ${prefix}.fastp.html --json ${prefix}.fastp.json
    """
  }
  // single end not yet supported
}