#!/usr/bin/env Rscript

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressMessages(library(ShortRead))


## expected barcodes
if(length(commandArgs(trailingOnly=TRUE))==3) {
  a <- commandArgs(trailingOnly=TRUE)[3]
} else {
  a <- "1-8"
}

a.split <- lapply(strsplit(strsplit(a, ",")[[1]], "-"), as.numeric)
a.split[lengths(a.split)==2] <- lapply(a.split[lengths(a.split)==2], function(x) seq(x[1], x[2]))
# a.split

## known barcode sequences
b <- c(AD001="ATCACG", AD002="CGATGT", AD004="TGACCA", AD006="GCCAAT", AD007="CAGATC", AD008="ACTTGA", AD010="TAGCTT", AD012="CTTGTA")

# only expected barcodes
bc0 <- b[unlist(a.split)]
# bc0

message("looking for ", length(bc0), " barcodes:")
for(i in names(bc0)) message(" - ", i, ": ", bc0[i])

# generate all barcode sequences with 1mm
per <- function(x) unlist(lapply(lapply(seq(nchar(x)), function(i) paste0(substr(x, 1, i-1), c("A", "C", "G", "T"), substr(x, i+1, 100))), setdiff, x))
bc0.per <- lapply(bc0, per)

bc1 <- unlist(bc0.per)
names(bc1) <- paste0(rep(names(bc0.per), lengths(bc0.per)), "-1mm")
# bc1

# merge all barcodes
bc <- c(bc0, bc1)
# bc


## read and write fastq files at the same time
f <- commandArgs(trailingOnly=TRUE)[1]
f.split <- strsplit(f, ",")[[1]]

# output file base
o <- commandArgs(trailingOnly=TRUE)[2]
#o <- gsub("_L001_R1_001.fastq.gz", "", f.split)

for(i in names(bc0)) {
  if(file.exists(paste0(o, "-", i, "_R1.fastq.gz")) | file.exists(paste0(o, "-", i, "_R2.fastq.gz"))) stop("Output file exists")
}

# matrix to collect stats per barcode
n <- matrix(0, length(bc0)+1, 2, dimnames = list(paste0(o, "-", c(names(bc0), "NA")), c("exact", "1mm")))

# loop over input files
for(f1 in f.split) {
  #message(f1)
  f2 <- sub("_R1", "_R2", f1)
  
  strm1 <- FastqStreamer(f1, n=1E6)  # 1M reads by default
  strm2 <- FastqStreamer(f2, n=1E6)
  
  message("processing file ", match(f1, f.split), "/", length(f.split), " (", f1, "): ", appendLF = FALSE)
  repeat {
    message(".", appendLF = FALSE)
    fq1 <- yield(strm1)
    if(length(fq1) == 0) break
    fq2 <- yield(strm2)
    
    # match to barcodes
    fq1.six <- as.vector(subseq(sread(fq1), 1, 6))
    fq1.bc <- names(bc)[match(fq1.six, bc)]
    
    # count barcode occurence
    n[, 1] <- n[, 1] + c(table(factor(fq1.bc, unique(names(bc))))[names(bc0)], sum(is.na(fq1.bc)))
    n[, 2] <- n[, 2] + c(table(factor(fq1.bc, unique(names(bc))))[paste0(names(bc0), "-1mm")], NA)
    
    # Remove cell barcodes (trim first 6 bases of R1). Added by Fabio 06/02/2022   
    #fq1 = narrow(fq1, start = 7)
    
    
    # write compressed fastq
    for(i in names(bc0)) {
      writeFastq(fq1[which(fq1.bc == i | fq1.bc == paste0(i, "-1mm"))], file = paste0(o, "-", i, "_R1.fastq.gz"), mode = "a", compress = TRUE)
      writeFastq(fq2[which(fq1.bc == i | fq1.bc == paste0(i, "-1mm"))], file = paste0(o, "-", i, "_R2.fastq.gz"), mode = "a", compress = TRUE)
    }
    writeFastq(fq1[which(is.na(fq1.bc))], file = paste0(o, "-", "NA", "_R1.fastq.gz"), mode = "a", compress = TRUE)
    writeFastq(fq2[which(is.na(fq1.bc))], file = paste0(o, "-", "NA", "_R2.fastq.gz"), mode = "a", compress = TRUE)
    
    invisible(gc())
  }
  message(" done")
  
  close(strm1)
  close(strm2)
}

# write summary
write.table(data.frame(n, percent=sprintf("%.1f%%", rowSums(n, na.rm = TRUE)/sum(n, na.rm = TRUE)*100)), file = paste0(o, ".bcsplit.txt"), sep="\t", quote=FALSE)
message(sum(n, na.rm = TRUE), " reads processed")