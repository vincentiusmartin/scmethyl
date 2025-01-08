
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCSPLIT } from '../modules/custom/bcsplit'
include { FASTP } from '../modules/fastp/fastp'
include { METHYLCTOOLS_FQCONV } from '../modules/methylctools/fqconv'
include { METHYLCTOOLS_BCALL } from '../modules/methylctools/bcall'
include { SAMBAMBA_MARKDUP } from '../modules/sambamba/markdup'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCMETHV1 {
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
    
  // Perform demultiplexing to per-cell processing
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
    
  FASTP(ch_demultiplexed) 
  FASTP.out.reads | METHYLCTOOLS_FQCONV | SAMBAMBA_MARKDUP | METHYLCTOOLS_BCALL
}
