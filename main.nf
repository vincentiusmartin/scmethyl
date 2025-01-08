#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SCMETHV1 } from './workflows/scmethv1'
include { SCMETHV2 } from './workflows/scmethv2'

workflow {
    if (params.version == "v2") {
        SCMETHV2()
    } else {
        SCMETHV1()
    }
}