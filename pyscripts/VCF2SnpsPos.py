#!/bin/env python3

import gzip
import argparse
import os

def get_arguments():
    """function to return the paths to the input files"""

    parser = argparse.ArgumentParser(description="Generate a table of the snps positions used in MatrixEQTL")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="VCF containing SNPs to be used in the MatrixEQTL analysis")

    args = parser.parse_args()
    vcf = args.vcf

    # make sure input files exist
    assert os.path.isfile(vcf), f"Could not locate {vcf}"

    return vcf

def make_output_name(vcf: str):
    """helper function to create output name"""

    bname = os.path.basename(vcf)

    if bname.endswith(".gz"):
        extension = ".vcf.gz"
    else:
        extension = ".vcf"

    outfile = bname.replace(extension, "_snpspos.tsv")

    return outfile

def parse_vcf(vcf: str):
    """Get the SNP positions and write them out to a tsv file"""
    
    # create output name
    outname = make_output_name(vcf)

    outfh = open(outname, 'w')
    outfh.write("snpid\tchr\tpos\n")

    # now for the filtering
    fh = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.strip('\n').split('\t')
        chrom = fields[0]
        snp_pos = int(fields[1])
        snp_id = f"{chrom}_{snp_pos}" # create a custom name for this SNP
        outline = f"{snp_id}\t{chrom}\t{snp_pos}\n"
        outfh.write(outline)

    # close IOs
    fh.close()
    outfh.close()
    

def main():
    """entry point to the pipeline"""

    # get the input
    vcf = get_arguments()

    # simply write to a table
    parse_vcf(vcf)

if __name__ == "__main__":
    main()