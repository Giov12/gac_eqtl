#!/bin/env python3

import gzip
import argparse
import os
import sys

def get_arguments():
    """function to return the paths to the htseq counts file and list of samples"""

    parser = argparse.ArgumentParser(description="Subsample a HTSeqs table for the samples provided in the samples file")
    parser.add_argument("-c", "--counts", required=True, type=str, help="htseq counts file")
    parser.add_argument("-s", "--samples", required=True, type=str, help="A single column list of the samples to select from")


    args = parser.parse_args()
    counts_file = args.counts
    samples_file = args.samples

    # make sure input file exist
    for f in [counts_file, samples_file]:
        assert os.path.isfile(f), f"Could not locate {f}"

    return counts_file, samples_file

def get_samples(samples_file: str):
    """get the samples from the samples file"""

    fh = gzip.open(samples_file, "rt") if samples_file.endswith(".gz") else open(samples_file, 'r')

    samples = set()

    for line in fh:
        if len(line) > 0:
            samples.add(line.strip())

    fh.close()

    return samples

def make_output_name(counts_file: str):
    """create the name of the output file"""

    bname = os.path.basename(counts_file)

    if ".tsv" in bname:
        extension = ".tsv"
    elif ".txt" in bname:
        extension = ".txt"
    else:
        print("Unsure of extension. Going to append only when making output file")
        extension = ""
    # if compressed
    if bname.endswith(".gz"):
        extension = extension + ".gz"

    outname = bname.replace(extension, "_subsampled.txt")

    return outname
    

def subsample_counts(counts_file: str, samples: set):
    """here is where we create the output (i.e., subsample counts)"""

    outname = make_output_name(counts_file)
    outfh = open(outname, 'w')

    fh = gzip.open(counts_file, "rt") if counts_file.endswith(".gz") else open(counts_file, 'r')

    # assume first line has the sample headers
    first_line = True

    for line in fh:
        fields = line.strip().split('\t')
        if first_line:
            sample_indices = []
            sample_columns = [fields[0]] # assume for gene IDs
            for i, sample in enumerate(fields):
                sample = sample.replace('D', '') # hard coding this in
                if sample in samples:
                    sample_indices.append(i)
                    sample_columns.append(sample)
            if len(sample_indices) == 0:
                msg = "Found 0 matching samples"
                sys.exit(msg)
            first_line = False
            outline = '\t'.join(sample_columns) + '\n'
            outfh.write(outline)
        else:
            outline = fields[0] # use gene name
            for i in sample_indices:
                outline += '\t' + fields[i]
            outfh.write(outline + '\n')

    fh.close()
    outfh.close()


def main():
    """entry point to the pipeline"""

    # get the inputs
    counts_file, samples_file = get_arguments()

    # get the samples
    samples = get_samples(samples_file)

    # subsample and write output
    subsample_counts(counts_file, samples)

if __name__ == "__main__":
    main()