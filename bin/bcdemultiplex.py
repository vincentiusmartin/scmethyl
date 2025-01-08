#!/usr/bin/env python

import fastools
import pandas as pd
import edlib
import argparse
import os
import gzip
from jit_open import Handle, Queue
from Bio import SeqIO
from gzip import open as gzip_open
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Demultiplexing single cell with barcode on R1")
    parser.add_argument('--fq1', type=str, nargs='+', 
                        help="List of forward FastQ file paths (R1 files)", required=True)
    parser.add_argument('--fq2', type=str, nargs='+', 
                        help="List of reverse FastQ file paths (R2 files)", required=True)
    parser.add_argument('--bc', type=str, 
                        help="Path to the barcode file", required=True)
    parser.add_argument('--sample_name', type=str, 
                        help="Name of the sample", required=True)
    parser.add_argument('--out_path', type=str, 
                        help="Output directory path (default: current directory)",
                        default=os.getcwd())
    parser.add_argument('--mismatch', type=int, 
                        help="Mismatch threshold for barcode match", default=1)
    args = parser.parse_args()
    return args

def open_files(path, samplename, barcodeids, queue):
    if not os.path.exists(path):
        os.path.mkdir(path)
    handles = {}
    for bcid in barcodeids:
        # hard coded for gz file, might need to be changed depending on needs
        handles['{}_R1'.format(bcid)] = Handle('{}/{}_{}_R1.fastq.gz'.format(path, samplename, bcid), queue, f_open=gzip_open)
        handles['{}_R2'.format(bcid)] = Handle('{}/{}_{}_R2.fastq.gz'.format(path, samplename, bcid), queue, f_open=gzip_open)
    return handles

def read_barcodes(file_path):
    with open(file_path, 'r') as f:
        sequences_dict = {line.strip().split()[1]: line.strip().split()[0] for line in f}
    return sequences_dict

def read_fasta(fastq_files):
    if isinstance(fastq_files, str):
        fastq_files = [fastq_files]
    for fastq_file in fastq_files:
        with gzip.open(fastq_file, "rt") as f:
            file_format = fastools.guess_file_format(f)
            reader = SeqIO.parse(f, file_format)
            for record in reader:
                yield record

def get_match(seq, barcodes, mismatch=1):
    for barcode in barcodes:
        alignment = edlib.align(seq, barcode)
        if alignment['editDistance'] <= mismatch:
            return(barcodes[barcode]),alignment['editDistance']
    return "NA",None

def main():
    # Parse and assign arguments to variables
    args = parse_args()
    fq1_paths = args.fq1
    fq2_paths = args.fq2
    bc_path = args.bc
    sample_name = args.sample_name
    out_path = args.out_path
    mismatch = args.mismatch

    # Read the barcodes
    barcodes = read_barcodes(bc_path)
    lenbc = len(list(barcodes.keys())[0])
    queue = Queue()
    handles = open_files(out_path, sample_name, list(barcodes.values())+["NA"], queue)

    # Process the fastq files
    with open(fq1_paths[0], "rb") as f:
        file_format = fastools.guess_file_format(f)

    fq1s = read_fasta(fq1_paths)
    fq2s = read_fasta(fq2_paths)

    counter = defaultdict(int)
    for record1, record2 in zip(fq1s, fq2s):
        cellbc = record1.seq[:lenbc]
        bcmatch,dist = get_match(cellbc, barcodes)
        counter[(bcmatch,dist)] += 1
        SeqIO.write(record1, handles[bcmatch+"_R1"], file_format)
        SeqIO.write(record2, handles[bcmatch+"_R2"], file_format)
    queue.flush()

    # Write the barcode counts to a csv file
    df = pd.DataFrame(counter.items(), columns=['idmm', 'count'])
    df[['id', 'mismatch']] = pd.DataFrame(df['idmm'].tolist(), index=df.index)
    df = df.drop(columns='idmm')
    df = df[['id', 'mismatch', 'count']]
    df['mismatch'] = df['mismatch'].astype('Int32')
    df = df.sort_values(by=['id','mismatch'])
    df['percentage'] = df['count'] / df['count'].sum() * 100
    df.to_csv(out_path + sample_name + "_bccounts.csv", index=False)

if __name__ == "__main__":
    main()
