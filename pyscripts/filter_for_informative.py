#!/bin/env python3

import gzip
import argparse
import os
import sys
from typing import TextIO
from collections import defaultdict

meta_file  = ''
meta_info  = dict()
vcf_file   = list()
eqtl_file  = ''
info_sites = dict() # key -> chr & int(bp), value Info Site Class
minSupport = 0.8

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

    global vcf_file, eqtl_file, meta_file, minSupport

    desc  = "Filter a table of eQTLs from stickleback for fixed and/or private common/white alleles"
    vhelp = "A vcf file with common and white stickleback"
    ehelp = "eQTL file with common and white stickleback"
    mhelp = "A file containing the population assignments for each sample"
    shelp = "Minimum fraction of samples needed per population to support the site"

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-v", "--vcf",     required=True,      type=str,   help=vhelp)
    parser.add_argument("-e", "--eqtls",   required=True,      type=str,   help=ehelp)
    parser.add_argument("-m", "--meta",    required=True,      type=str,   help=mhelp)
    parser.add_argument("-s", "--support", default=minSupport, type=float, help=shelp)

    args       = parser.parse_args()
    vcf_file   = args.vcf
    eqtl_file  = args.eqtls
    meta_file  = args.meta
    minSupport = args.support

    # make sure input files exist
    assert os.path.isfile(vcf_file),  f"Could not locate {vcf_file}"
    assert os.path.isfile(eqtl_file), f"Could not locate {eqtl_file}"
    assert os.path.isfile(meta_file), f"Could not locate {meta_file}"

    if (minSupport > 1.0 or minSupport < 0):
        msg = "Error: --support must be between 0 and 1.0"
        sys.exit(msg)

    return 0

def open_file(infile: str) -> TextIO:
    """"helper function to open a file"""
    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')

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

    fh = open_file(meta_file)

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

    global vcf_file, info_sites, meta_info, minSupport

    fh         = open_file(vcf_file)
    ecotypes   = list()
    samples    = list()
    pop_counts = defaultdict(int)

    #
    # ecotypes will be encoded as
    # Common -> c, Whites -> w
    #

    #
    # Hard code the sex chrom IDs since those
    # variants were only called using one of the
    # sexes (i.e., males all have ./. in X chrom &
    # females all have ./. on Y)
    # 
    XCHROM = "chrXIX"
    YCHROM = "chrY"

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
                    samples.append(sample)
                    if (sample not in meta_info):
                        # we only care about the samples provided in the
                        # meta data
                        ecotypes.append('S') # S == Skip
                        continue
                    ecotype = meta_info[sample][0]
                    ecotypes.append(ecotype)
                    pop_counts[ecotype] += 1
                if (len(ecotypes) != len(meta_info)):
                    msg = "samples in vcf and meta information are not equal. Unreliable"
                    sys.exit(msg)
            continue
        fields   = line.strip().split('\t')
        chrom    = fields[0]
        wAlleles = set()
        cAlleles = set()
        counts   = {p : cnt for p, cnt in pop_counts.items()}

        # create an allele mapping
        nucs = [fields[3]]
        for nuc in fields[4].split(','):
            nucs.append(nuc)

        for i in range(9, len(fields)):
            ecotype = ecotypes[i - 9]
            sample  = samples[i - 9] 
            if (ecotype == 'S'):
                continue # skip
            info    = fields[i]
            idx     = info.find(':')
            geno    = info[:idx]
            alleles = geno.split('/')
            if (len(alleles) != 2 and geno.count('|') > 0):
                alleles = geno.split('|')
            haploid = False
            if (len(alleles) != 2):
                if (chrom != YCHROM and chrom != XCHROM):
                    msg = f"Encountered an malformated entry in the vcf\n" + \
                          f"Offending line: {line}"
                    sys.exit(msg)
                else:
                    haploid = True
            if (alleles[0] == '.'):
                # avoid penalizing sites where sex chroms were only called
                # for a specific sex
                if (chrom != YCHROM and chrom != XCHROM):
                    counts[ecotype] -= 1
                elif (chrom == YCHROM and sample.startswith('M')):
                    counts[ecotype] -= 1
                elif (chrom == XCHROM and sample.startswith('F')):
                    counts[ecotype] -= 1
                continue # skip missing calls
            if (ecotype == 'c'):
                if (haploid and len(alleles) == 1):
                    cAlleles.add(alleles[0])
                else:
                    cAlleles.add(alleles[0])
                    cAlleles.add(alleles[1])
            else:
                if (haploid and len(alleles) == 1):
                    wAlleles.add(alleles[0])
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
                # now test for min support
                for p, cnt in counts.items():
                    tot = pop_counts[p]
                    if (cnt / tot < minSupport):
                        valid = False
                        break
            if (valid):
                chrom  = fields[0].replace("chr", '') if fields[0].startswith("chr") else fields[0]
                offset = 71 if chrom == "XVIII" else 0
                bp     = int(fields[1]) - offset
                if (bp > 0): # this bp position exists in the Ensembl assembly
                    site = InfoSite()
                    site.common_alleles = set([nucs[int(n)] for n in cAlleles])
                    site.white_alleles  = set([nucs[int(n)] for n in wAlleles])
                    info_sites[(chrom, bp)] = site

    fh.close()

    if (len(info_sites) == 0):
        msg = f"Did not find any informative sites between commons and whites"
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
    """check the eQTL and see if we can find any sites that overlap with the informative sites"""

    global info_sites, eqtl_file

    fh  = open_file(eqtl_file)
    ofh = open(make_output_name(), 'w')
    cnt = 0 # count of number of records written
    lnm = 0 # line number

    for line in fh:
        if (len(line) == 0):
            continue
        line   = line.rstrip("\r\n\t ")
        fields = line.split('\t')
        last   = fields[-1] + '\n'
        lnm   += 1
        if (lnm == 1): # should be the header
            # add a new column
            fields[-1] = "Allele.Info"
            fields.append(last)
            line = '\t'.join(fields)
            ofh.write(line)
        else:
            chrom  = fields[1]
            bp     = int(fields[2])
            if ((chrom, bp) in info_sites):
                site = info_sites[(chrom, bp)]
                # add a new column
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

    # load the meta information for each sample
    load_population_assignments()

    # find the informative sites
    get_informative_sites()

    # now filter the eqtls
    filter_eqtls()

if __name__ == "__main__":
    main()