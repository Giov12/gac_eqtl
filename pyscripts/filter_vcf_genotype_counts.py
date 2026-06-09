#!/bin/env python3

import gzip
import argparse
import os
from   collections import defaultdict

vcf_file = ''
min_idv  = 0

def get_arguments() -> int:
    """function to set the global variables"""

    global vcf_file, min_idv

    parser = argparse.ArgumentParser(description="Filter VCF based on the number of minimum individuals per genotype")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="VCF file")
    parser.add_argument("-c", "--count", required=True, type=int, help="Minimum number of individuals with a genotype")

    args     = parser.parse_args()
    vcf_file = args.vcf
    min_idv  = args.count

    # make sure input file exist
    assert os.path.isfile(vcf_file), f"Could not locate {vcf_file}"

    # make sure max distance is legal
    assert min_idv >= 0, "--count cannot be less than 0"

    return 0

def filter_vcf():
    """Work horse function of this program"""

    global vcf_file, min_idv

    # create output name
    replacement = ".vcf.gz" if vcf_file.endswith(".gz") else ".vcf"
    outname     = "Filtered." + os.path.basename(vcf_file).replace(replacement, '') + f".{min_idv}genotypes.vcf"
    outfh       = open(outname, 'w')
    fh          = gzip.open(vcf_file, "rt") if vcf_file.endswith(".gz") else open(vcf_file, 'r')
    kept        = 0
    filtered    = 0

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            outfh.write(line)
            continue

        else:
            fields = line.strip('\n').split('\t')
            counts = defaultdict(int)
            for i in range(9, len(fields)):
                info = fields[i]
                j    = info.find(':')
                gt   = info[:j]
                if (gt != "./."):
                    counts[gt] += 1
            keep = True
            if (len(counts) > 0):
                for count in counts.values():
                    if (count < min_idv):
                        keep = False
                        break
            else:
                keep = False
            if (keep):
                outfh.write(line)
                kept += 1
            else:
                filtered += 1

    # close IOs
    fh.close()
    outfh.close()

    print(f"Filtered {filtered} variant positions. Retained {kept} variants")
    

def main():
    """entry point to the pipeline"""

    # get the inputs
    get_arguments()

    # now to create the new VCF
    filter_vcf()

if __name__ == "__main__":
    main()