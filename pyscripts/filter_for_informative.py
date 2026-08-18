#!/bin/env python3

import gzip
import argparse
import os
import sys

meta_file  = ''
meta_info  = dict()
vcf_file   = list()
eqtl_file  = ''
info_sites = dict() # key -> chr & int(bp), value Info Site Class

class InfoSite:
    def __init__(self) -> None:
        self.common_alleles = set()
        self.white_alleles  = set()

    def get_repr(self) -> str:
        str1 = "Common (" + ','.join(self.common_alleles) + ')'
        str2 = "Whites (" + ','.join(self.white_alleles) + ')'
        return str1 + ", " + str2
   
def get_arguments() -> int:
    """function to set the global variables"""

    global vcf_file, eqtl_file, meta_file

    parser = argparse.ArgumentParser(description="Filter a table of eQTLs from stickleback for fixed and/or private common/white alleles")
    parser.add_argument("-v", "--vcf",   required=True, type=str, help="a vcf files with common and white stickleback")
    parser.add_argument("-e", "--eqtls", required=True, type=str, help="eqtl file with common and white stickleback")
    parser.add_argument("-m", "--meta",  required=True, type=str, help="a file containing the population assignments for each sample")

    args      = parser.parse_args()
    vcf_file  = args.vcf
    eqtl_file = args.eqtls
    meta_file = args.meta

    # make sure input files exist
    assert os.path.isfile(vcf_file), f"Could not locate {vcf_file}"
    assert os.path.isfile(eqtl_file), f"Could not locate {eqtl_file}"
    assert os.path.isfile(meta_file), f"Could not locate {meta_file}"

    return 0

def load_population_assignments() -> int:
    """load in the ecotype assignments for each of the samples"""

    global meta_file, meta_info

    #
    # expected file structure
    # sample \t ecotype
    #
    # for this analysis, we will
    # create a mapping of sample -> ecotype
    #

    fh = gzip.open(meta_file, "rt") if meta_file.endswith(".gz") else open(meta_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 2, f"Expected 2 columns. Look at line {line}"
        sample  = fields[0]
        ecotype = fields[1].strip()
        meta_info[sample] = ecotype

    fh.close()

    print(f"Loaded {len(meta_info)} samples from {meta_file}")

    return 0

def get_informative_sites() -> int:
    """Find the load the informative sites between commons vs whites"""

    global vcf_file, info_sites, meta_info

    fh          = gzip.open(vcf_file, "rt") if vcf_file.endswith(".gz") else open(vcf_file, 'r')
    phenotypes  = list()

    #
    # ecotypes will be encoded as
    # Common -> c, Whites -> w
    #

    #  
    # !!!!!!!!!!!!!
    # NOTE: chromosome XVIII has a stretch of 70 N's
    # at the very beginning, so we need to offset
    # SNPs on this chromosome in order to match
    # the SNPs on Ensembl v113
    #
    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields = line.strip().split('\t')
                for i in range(9, len(fields)):
                    sample = fields[i]
                    if (sample not in meta_info):
                        # we only care about the samples provided in the
                        # meta data
                        phenotypes.append('S') # S == Skip
                        continue
                    ecotype = meta_info[sample][0]
                    phenotypes.append(ecotype)
                if (len(phenotypes) != len(meta_file)):
                    msg = "samples in vcf and meta information are not equal. Unreliable"
                    sys.exit(msg)
            continue
        fields   = line.strip().split('\t')
        wAlleles = set()
        cAlleles = set()

        for i in range(9, len(fields)):
            pheno   = phenotypes[i - 9]
            if (pheno == 'S'):
                continue # skip
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
                chrom  = fields[0].replace("chr", '') if fields[0].startswith("chr") else fields[0]
                offset = 71 if chrom == "XVIII" else 0
                bp     = int(fields[1]) - offset
                if (bp > 0): # this bp position exists in the Ensembl assembly
                    site = InfoSite()
                    site.common_alleles = cAlleles
                    site.white_alleles  = wAlleles
                    info_sites[(chrom, bp)] = site

    fh.close()

    if (len(info_sites) == 0):
        msg = f"Did not find any fixed differences between commons and whites"
        sys.exit(msg)
    
    print(f"Found a total of {len(info_sites)} informative sites between commons and whites")

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

    global info_sites, eqtl_file

    fh  = gzip.open(eqtl_file, "rt") if eqtl_file.endswith(".gz") else open(eqtl_file, 'r')
    ofh = open(make_output_name(), 'w')
    cnt = 0 # count of number of records written
    lnm = 0 # line number

    for line in fh:
        if (len(line) == 0):
            continue
        fields = line.split('\t')
        lnm += 1
        if (lnm == 1): # should be the header
            # add a new column
            last  = fields[-1]
            label = "Allele.Info"
            fields[-1] = label
            fields.append(last)
            line = '\t'.join(fields)
            ofh.write(line)
        else:
            chrom  = fields[0]
            bp     = int(fields[1])
            if ((chrom, bp) in info_sites):
                site = info_sites[(chrom, bp)]
                # add a new column
                last = fields[-1]
                # add the new column
                fields[-1] = site.get_repr()
                fields.append(last) # re-add this column
                line = '\t'.join(fields)
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

    # find the informative sites
    get_informative_sites()

    # now filter the eqtls
    filter_eqtls()

if __name__ == "__main__":
    main()