
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { BCDEMUL } from '../modules/custom/bcdemultiplex'
include { BCSPLIT } from '../modules/custom/bcsplit'
include { FASTP } from '../modules/fastp/fastp'
include { BISCUIT_ALIGN } from '../modules/biscuit/align'
include { SAMBAMBA_MARKDUP } from '../modules/sambamba/markdup'
include { METHYLDACKEL } from '../modules/methyldackel/methyldackel'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCMETHV2 {
  Channel
    .fromPath(params.input)
    .splitCsv(skip: 1)
    .map {
        meta, fastq_1, fastq_2 ->
        if (!fastq_2) {
            return [ [id:meta, single_end:true ], [ fastq_1 ] ]
        } else {
            return [ [id:meta, single_end:false ], [ fastq_1, fastq_2 ] ]
        }
    }
    .groupTuple(by: [0])
    .map { meta, fastq -> [meta,fastq.flatten()] }
    .set { ch_fastq }
    
    BCSPLIT(ch_fastq)
      .reads
      .transpose()
      .map{
        meta,fastq -> 
          def newid = (fastq.getSimpleName() =~ /(.*)_R[12]/)[0][1]
          return [[id:newid, single_end:meta.single_end], [fastq]]
        }
      .groupTuple(by: [0])
      .map { meta, fastq -> [meta,fastq.flatten()] }
      .set{ ch_demultiplexed }

    // MODULE: BBDUK_TRIM
    FASTP(ch_demultiplexed)
    // need to create biscuit index before running
    FASTP.out.reads | BISCUIT_ALIGN | SAMBAMBA_MARKDUP
    METHYLDACKEL(SAMBAMBA_MARKDUP.out.reads)
}

