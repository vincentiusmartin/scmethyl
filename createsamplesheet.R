library(tidyverse)
library(data.table)

datadir <- "/research/groups/northcgrp/projects/northcgrp_hartwell/common/illumina"
metadatapath <- "/research/groups/northcgrp/home/common/Vincentius/snmcseq/scWGBS_results/CellAnnotations_v0.txt"
# in general, we always need directory name and sample name
dircol <- "Sequencing_Round"
samplecol <- "Multiplexed_file"

search_files <- function(indir, samples, pattern) {
  indir <- indir[[1]]
  patterns <- sprintf(pattern,samples)
  allfiles <- list.files(indir, recursive=TRUE, full.names=TRUE)
  matches <- lapply(patterns, function(p) allfiles[grepl(p,allfiles)])
  return(matches)
}

metadata <- fread(metadatapath)
samplesheet <- distinct(metadata,!!sym(dircol),!!sym(samplecol)) %>%
  rename(seqdir=!!sym(dircol),sample=!!sym(samplecol)) %>%
  mutate(seqdir = file.path(datadir,seqdir)) %>%
  group_by(seqdir) %>%
  mutate(fastq_1 = search_files(seqdir, sample, "*%s(.*)_R1_*")) %>%
  ungroup() %>%
  unnest(fastq_1) %>%
  select(-seqdir) %>%
  distinct() %>%
  mutate(fastq_2 = str_replace(fastq_1,"_R1_", "_R2_"),
         fastq_2 = if_else(file.exists(fastq_2),fastq_2,NA)) 
fwrite(samplesheet,"samplesheet.csv")
