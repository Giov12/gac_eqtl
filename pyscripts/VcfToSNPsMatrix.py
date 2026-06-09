#!/bin/env python3

import gzip
import argparse
import os

def get_arguments():
    """function to return the paths to the vcf file to reformat"""

    parser = argparse.ArgumentParser(description="Reformat a VCF file to match a SNP matrix for MatrixEQTL")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="VCF file")


    args = parser.parse_args()
    vcf = args.vcf

    # make sure input file exist
    assert os.path.isfile(vcf), f"Could not location {vcf}"

    return vcf


def reformat_vcf(vcf: str):
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
            samples = fields[9:]
            header_line = "snpid\t"
            header_line = header_line + '\t'.join(samples) + '\n'
            outfh.write(header_line)
        else:
            fields = line.strip('\n').split('\t')
            genos = [] # will be used to create the outline
            snp_id = fields[0] + '_' + fields[1]
            genos.append(snp_id)
            genotypes = fields[9:]
            for geno in genotypes:
                geno_info = geno.split(':')
                genos.append(geno_key[geno_info[0]])
            outline = '\t'.join(genos) + '\n'
            outfh.write(outline)

    # close IOs
    fh.close()
    outfh.close()
    

def main():
    """entry point to the pipeline"""

    # get the input
    vcf = get_arguments()

    # reformat
    reformat_vcf(vcf)

    # creatte output
if __name__ == "__main__":
    main()