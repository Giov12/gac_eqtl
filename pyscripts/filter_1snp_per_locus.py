#!/bin/env python3

import gzip
import argparse
import os

vcf_file = ''

class Locus:
    def __init__(self, id_: str) -> None:
        self.id      = id_
        self.records = list()

    def size(self) -> int:
        return len(self.records)
    
    def add_record(self, record: str) -> int:
        self.records.append(record)
        return 0
    
    def get_1_snp(self) -> str:

        #
        # no need to do a bunch of work
        # if there's only 1 option
        #
        if (len(self.records) == 1):
            return self.records[0]

        counts = list()
        
        #
        # get the record with the highest NS
        # if tied, choose the first one
        #
        cur_max = -1
        for record in self.records:
            fields = record.split('\t')
            info   = fields[7]
            idx    = 3
            jdx    = -1
            for i in range(len(info)):
                if (info[i] == ';'):
                    jdx = i
                    break
            count = int(info[idx:jdx])
            if (count > cur_max):
                cur_max = count
            counts.append(count)

        for i in range(len(counts)):
            count = counts[i]
            if (count == cur_max):
                return self.records[i]

def get_arguments() -> int:
    """function to set the global variables"""

    global vcf_file

    parser = argparse.ArgumentParser(description="Filter a vcf produced by populations to keep 1 SNP per locus")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="populations.snps.vcf.gz file")

    args     = parser.parse_args()
    vcf_file = args.vcf

    # make sure input file exist
    assert os.path.isfile(vcf_file), f"Could not locate {vcf_file}"

    return 0

def filter_vcf():
    """Work horse function of this program"""

    global vcf_file

    # create output name
    gzipped  = vcf_file.endswith(".gz")
    outname  = "1SNP.per.locus." + os.path.basename(vcf_file) 
    outfh    = gzip.open(outname, "wt") if gzipped else open(outname, 'w')
    fh       = gzip.open(vcf_file, "rt") if gzipped else open(vcf_file, 'r')
    curLocus = None
    kept     = 0
    seen     = 0

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            outfh.write(line)
            continue

        else:
            seen   += 1
            fields  = line.split('\t')
            snpID   = fields[2]
            idx     = snpID.find(':')
            locusID = snpID[:idx]
            if (curLocus == None):
                curLocus = Locus(locusID)
            if (curLocus.id == locusID):
                curLocus.add_record(line)
            else:
                outline = curLocus.get_1_snp()
                outfh.write(outline)
                curLocus = Locus(locusID)
                curLocus.add_record(line)
                kept    += 1

    fh.close()

    # get last locus
    outline = curLocus.get_1_snp()
    outfh.write(outline)
    outfh.close()

    print(f"Processed {seen} variant positions. Retained {kept} variants")
    

def main():
    """entry point to the pipeline"""

    # get the inputs
    get_arguments()

    # now to create the new VCF
    filter_vcf()

if __name__ == "__main__":
    main()