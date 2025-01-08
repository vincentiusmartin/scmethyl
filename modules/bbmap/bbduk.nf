process BBDUK_TRIM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/bbduk", mode: 'copy'

    module "bbmap"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // set number of cores
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }

    def args_list = ["ktrim=r", "k=${params.bbduk_k}", "hdist=${params.hdist}", 
                    "mink=${params.mink}", "ref=${params.adapter_ref}", "rcomp=f",
                    "ordered=t","tbo=t","tpe=t", "threads=$cores"]

    if (meta.single_end) {
        """
        bbduk.sh \\
            ${args_list.join(' ')} \\
            in=$reads \\
            out=${prefix}.trimmed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        END_VERSIONS
        """
    
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

        bbduk.sh \\
            ${args_list.join(' ')} \\
            in1=${prefix}_1.fastq.gz in2=${prefix}_2.fastq.gz \\
            out1=${prefix}_1.trimmed.fastq.gz out2=${prefix}_2.trimmed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        END_VERSIONS
        """
    }

}