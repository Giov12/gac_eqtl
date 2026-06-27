#!/bin/env python3

import gzip
import argparse
import os
import sys

vcf_files     = list()
eqtl_file     = ''
info_sites    = set() # tuples of chr & bp
   
def get_arguments() -> int:
    """function to set the global variables"""

    global vcf_files, eqtl_file

    parser = argparse.ArgumentParser(description="Filter a table of eQTLs from stickleback for fixed common/white alleles")
    parser.add_argument("-v", "--vcfs",   required=True, nargs='*', help="vcf files with common and white stickleback")
    parser.add_argument("-e", "--eqtls", required=True, type=str, help="eqtl file with common and white stickleback")

    args      = parser.parse_args()
    vcfs      = args.vcfs
    eqtl_file = args.eqtls

    # make sure input files exist
    assert os.path.isfile(eqtl_file), f"Could not locate {eqtl_file}"
    assert len(vcfs) > 0, "At least one vcf file must be provided"

    for vcf in vcfs:
        assert os.path.isfile(vcf), f"Could not locate {vcf}"
        vcf_files.append(vcf)

    return 0

def get_fixed_sites() -> int:
    """Find the load the informative sites between commons vs whites"""

    global vcf_files, info_sites

    #
    # Common samples have a prefix of CB
    # White samples have a prefix of CL
    #

    fh         = gzip.open(vcf_file, "rt") if vcf_file.endswith(".gz") else open(vcf_file, 'r')
    phenotypes = list() # w -> whites, c -> commons

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields = line.strip().split('\t')
                cnt    = 0
                for i in range(9, len(fields)):
                    cnt += 1
                    if (fields[i].startswith("CB")):
                        phenotypes.append('c')
                    elif (fields[i].startswith("CL")):
                        phenotypes.append('w')
                if (len(phenotypes) != cnt):
                    msg = "vcf contained samples without the CB/CL prefix. Unreliable"
                    sys.exit(msg)
            continue
        fields   = line.strip().split('\t')
        wAlleles = set()
        cAlleles = set()

        for i in range(9, len(fields)):
            pheno   = phenotypes[i - 9]
            info    = fields[i]
            idx     = info.find(':')
            geno    = info[:idx]
            alleles = geno.split('/')
            if (len(alleles) != 2 and geno.count('|') > 0):
                alleles = geno.split('|')
            if (len(alleles) != 2):
                msg = f"Encountered an ilformated entry in the vcf\n" + \
                      f"Offending line: {line}"
                sys.exit(msg)
            if (alleles[0] == '.'):
                continue # skip missing calls
            if (pheno == 'c'):
                cAlleles.add(alleles[0])
                cAlleles.add(alleles[1])
            else:
                wAlleles.add(alleles[0])
                wAlleles.add(alleles[1])

        # now check that they don't match
        if (len(wAlleles) > 0 and len(cAlleles) > 0):
            # if no alleles are shared
            valid = wAlleles.isdisjoint(cAlleles)
            if (valid == False): # else, check for private alleles
                if (len(wAlleles.difference(cAlleles)) > 0 or len(cAlleles.difference(wAlleles)) > 0):
                    valid = True
            if (valid):
                chrom = fields[0].replace("chr", '') if fields[0].startswith("chr") else fields[0]
                bp    = fields[1]
                fixed_alleles.add((chrom, bp))

    fh.close()

    if (len(fixed_alleles) == 0):
        msg = f"Did not find any fixed differences between commons and whites"
        sys.exit(msg)
    
    print(f"Found a total of {len(fixed_alleles)} fixed sites between commons and whites")

    return 0

def make_output_name() -> str:
    """create an output name for the eQTL dataset"""

    global eqtl_file

    bname = os.path.basename(eqtl_file)
    parts = list()

    for part in bname.split('.'):
        if (part == "tsv"):
            parts.append("filtered")
        parts.append(part)

    return '.'.join(parts)

    
def filter_eqtls() -> int:
    """check the eQTL and see if it is a fixed difference"""

    global fixed_alleles, eqtl_file

    fh  = gzip.open(eqtl_file, "rt") if eqtl_file.endswith(".gz") else open(eqtl_file, 'r')
    ofh = open(make_output_name(), 'w')
    cnt = 0 # count of number of records written
    lnm = 0 # line number

    for line in fh:
        if (len(line) == 0):
            continue
        lnm += 1
        if (lnm == 1): # should be the header
            ofh.write(line)
        else:
            fields = line.split('\t')
            chrom  = fields[0]
            bp     = fields[1]
            if ((chrom, bp) in fixed_alleles):
                ofh.write(line)
                cnt += 1

    fh.close()
    ofh.close()

    print(f"Encountered {lnm - 1} eQTL records. Kept {cnt} records")

    return 0

def main():
    """entry point to the pipeline"""

    # get the inputs
    get_arguments()

    # get the fixed differences
    get_fixed_sites()

    # now filter the eqtls
    filter_eqtls()

if __name__ == "__main__":
    main()