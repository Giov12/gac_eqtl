#!/bin/env python3

import gzip
import argparse
import os

def get_arguments():
    """function to return the paths to the vcf file to reformat and the file with the names of the samples to select"""

    parser = argparse.ArgumentParser(description="Reformat a VCF file to match a SNP matrix for MatrixEQTL and only select the samples provided")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="VCF file")
    parser.add_argument("-s", "--samples", required=True, type=str, help="A single column list of the samples to select from")


    args = parser.parse_args()
    vcf = args.vcf
    samples_file = args.samples

    # make sure input file exist
    for f in [vcf, samples_file]:
        assert os.path.isfile(f), f"Could not location {f}"

    return vcf, samples_file

def get_samples(samples_file: str):
    """get the samples from the samples file"""

    fh = gzip.open(samples_file, "rt") if samples_file.endswith(".gz") else open(samples_file, 'r')

    samples = set()

    for line in fh:
        if len(line) > 0:
            samples.add(line.strip())

    fh.close()

    return samples


def reformat_vcf(vcf: str, samples: set):
    """write out a matrix file to load into MatrixEQTL"""

    """
    0 = Homo Ref | 0/0
    1 = Hetero | 0/1 | 1/0
    2 = Homo Alt | 1/1
    """
    
    # create map for easy conversion
    geno_key = {"0/0" : '0', "0/1" : '1', "1/0" : '1', "1/1" : '2', './.' : 'NA'}

    # open & get extension
    if vcf.endswith(".gz"):
        fh = gzip.open(vcf, "rt")
        extension = ".vcf.gz"
    else:
        fh = open(vcf, 'r')
        extension = ".vcf"

    outname = os.path.basename(vcf).replace(extension, "_SNPsMatrix.tsv")
    outfh = open(outname, 'w')

    # now to parse through
    for line in fh:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM\t"):
            fields = line.strip('\n').split('\t')
            samples_key = []
            samples_columns = []
            header_line = "snpid\t"
            for i, col in enumerate(fields):
                if col in samples:
                    samples_key.append(i)
                    samples_columns.append(col)
            header_line = header_line + '\t'.join(samples_columns) + '\n'
            outfh.write(header_line)
        else:
            fields = line.strip('\n').split('\t')
            genos = [] # will be used to create the outline
            snp_id = fields[0] + '_' + fields[1]
            genos.append(snp_id)
            for sample_col in samples_key:
                genotype = fields[sample_col]
                geno_info = genotype.split(':')
                genos.append(geno_key[geno_info[0]])
            outline = '\t'.join(genos) + '\n'
            outfh.write(outline)

    # close IOs
    fh.close()
    outfh.close()
    

def main():
    """entry point to the pipeline"""

    # get the inputs
    vcf, samples_file = get_arguments()

    # get the samples
    samples = get_samples(samples_file)

    # subsample and reformat
    reformat_vcf(vcf, samples)

    # creatte output
if __name__ == "__main__":
    main()