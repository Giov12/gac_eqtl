#!/bin/env python3

import gzip
import argparse
import os
import sys

targets  = set()
geneFile = ''
dataFile = ''

def get_arguments() -> int:
    """function to return the paths to the target genes and the combined Counts/Genotypes tsv file"""

    global geneFile, dataFile

    parser = argparse.ArgumentParser(description="Get the records for target genes in the full expression/genotypes file")
    parser.add_argument("-g", "--genes", required=True, type=str, help="list of genes to get")
    parser.add_argument("-d", "--data", required=True, type=str, help="Full expression/genotypes file")

    args     = parser.parse_args()
    geneFile = args.genes
    dataFile = args.data

    # make sure input file exist
    for f in [geneFile, dataFile]:
        assert os.path.isfile(f), f"Could not locate {f}"

    return 0

def load_targets() -> int:
    """load the target genes"""

    global targets, geneFile

    fh = gzip.open(geneFile, "rt") if geneFile.endswith(".gz") else open(geneFile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        target = line.strip()
        targets.add(target)

    fh.close()

    return 0

def get_targets() -> int:
    """parse the tsv file for the target genes"""

    global targets, dataFile

    fh  = gzip.open(dataFile, "rt") if dataFile.endswith(".gz") else open(dataFile, 'r')
    ofh = open("targets.tsv", 'w')
    hdr = False # header has been written

    for line in fh:
        if (len(line) == 0):
            continue
        if (hdr == False):
            ofh.write(line)
            hdr = True
            continue
        fields = line.split('\t')
        gene   = fields[6]
        if (gene in targets):
            ofh.write(line)

    fh.close()
    ofh.close()

    return 0

def main():
    """entry point to the pipeline"""

    # get the inputs
    get_arguments()

    # load the target genes
    load_targets()

    # now parse
    get_targets()

if __name__ == "__main__":
    main()