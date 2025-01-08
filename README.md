# scmethyl

Single cell methylation pipeline from fastq to methylation call.

## Running pipeline

`nextflow run main.nf --input </path/to/samplesheet.csv> --version <v1/v2>`

Versions:

1. v1: using methylCtools
2. v2: using MethylDackel

Input samplesheet is a path to a csv file containing following fields:

| sample | fastq1 | fastq2 |
| --- | --- | --- |
| sample1 | /path/to/sample1_L001_R1.fastq.gz | /path/to/sample1_L001_R2.fastq.gz |
| sample1 | /path/to/sample1_L002_R1.fastq.gz | /path/to/sample1_L002_R2.fastq.gz |
| sample2 | /path/to/sample2_L001_R1.fastq.gz | /path/to/sample2_L001_R2.fastq.gz |
| sample2 | /path/to/sample2_L002_R1.fastq.gz | /path/to/sample2_L002_R2.fastq.gz |
