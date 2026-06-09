#!/bin/env python3

import gzip
import argparse
import os
from   collections import defaultdict

#
# globbal members
#
tgenes = set()
tsites = dict() # will be a str -> dict{str -> list[str]} | chr -> basepair : list[gene_ids]
gexprs = dict() # will be str -> dict{str : float} | individual -> gene : gene expression 


def get_arguments() -> tuple[str, str, str]:
    """function to return the paths to the input files"""

    parser = argparse.ArgumentParser(description="Generate a table of genotypes and target gene expression")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="VCF with only the target individuals")
    parser.add_argument("-t", "--targets", required=True, type=str, help="A .tsv file containing the chr and basepair position of the target sites along with the target gene_ids")
    parser.add_argument("-c", "--counts", required=True, type=str, help="Normalized read counts .tsv file")

    args = parser.parse_args()
    vcf = args.vcf
    tgt = args.targets
    cnt = args.counts


    # make sure input files exist
    assert os.path.isfile(vcf), f"Could not locate {vcf}"
    assert os.path.isfile(tgt), f"Could not locate {tgt}"
    assert os.path.isfile(cnt), f"Could not locate {cnt}"

    return (vcf, tgt, cnt)

def get_targets(tgt: str) -> int:
    """function to add the targets to our global members"""

    global tgenes
    global tsites

    fh = gzip.open(tgt, "rt") if tgt.endswith(".gz") else open(tgt, 'r')

    #
    # assuming the following format
    #   Type<tab>CHR<tab>BP<eGENE><tab>pvalue<tab>BH
    #
    ln = 0 
    for line in fh:
        ln += 1
        
        if (len(line) == 0 or line[0] == '#'):
            continue
        
        fields = line.strip().split('\t')
        assert len(fields) == 6, f"Encountered an unexpecting number at columns at line {ln}\n" + \
                                 f"Expected 6, found {len(fields)}"
        
        ety = fields[0] # eSNP type (cis || trans)
        chr = fields[1]
        bp  = fields[2]
        ge  = fields[3]
        pva = fields[4]
        bh  = fields[5]

        tgenes.add(ge) # add the gene
        if (chr not in tsites):
            tsites[chr] = defaultdict(list)
        tsites[chr][bp].append((ge, ety, pva, bh))

    fh.close()

    return 0

def get_expression(cnt: str) -> int:
    """get the expression of each target gene for each individual"""

    global gexprs
    global tgenes

    fh      = gzip.open(cnt, "rt") if cnt.endswith(".gz") else open(cnt, 'r')
    samples = []

    for line in fh:
        if (len(line) == 0 or (line[0] == '#')):
            continue
        fields = line.strip().split('\t')
        n      = len(fields)
        if (fields[0] == "gene_id"):
            for i in range(1, n):
                gexprs[fields[i]] = dict() # will hold gene expressions for target genes
                samples.append(fields[i])
            continue
        ge = fields[0]
        if (ge not in tgenes):
            continue
        for i in range(1, n):
            gexprs[samples[i - 1]][ge] = float(fields[i])

    fh.close()

    tgenes.clear() # we no longer need this in memory

    return 0

def parse_vcf(vcf: str) -> int:
    """core function to write out the genotype/expression data per individual"""
    
    global gexprs
    global tsites

    # create output name
    outname = "Deg.eQTL.GenoExpr.tsv"

    outfh   = open(outname, 'w')
    header  = "Chr\tBP\tType\tpvalue\tBH\tgene_id\t"  # will be written out after getting sample order
    samples = []

    # now for the filtering
    fh = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')

    for line in fh:
        if (len(line) == 0):
            continue
        if (line[0] == '#'):
            if (line.startswith("#CHROM")):
                fields  = line.strip().split('\t')
                samples = fields[9:]
                header  = header + '\t'.join(samples) + "\tGeno.Expr.Avgs"
                outfh.write(header)
                outfh.write('\n')
            continue
        fields = line.strip('\n').split('\t')
        chrom  = fields[0]
        bp     = fields[1]
        refGT  = fields[3] # reference allele
        altGt  = fields[4] # alternative allele
        genos  = [refGT, altGt]

        if ((chrom in tsites) and (bp in tsites[chrom])):
            genotypes = []

            for i in range(9, len(fields)):
                subfields = fields[i].split(':')
                genotype  = subfields[0]
                allele1   = genotype[0]
                allele2   = genotype[-1]
                # get the genotype
                if (allele1 != '.'):
                    if (allele1 == '0'):
                        allele1 = genos[0]
                    else:
                        allele1 = genos[1]
                if (allele2 != '.'):
                    if (allele2 == '0'):
                        allele2 = genos[0]
                    else:
                        allele2 = genos[1]
                genotype = f"{allele1}/{allele2}"
                genotypes.append(genotype)
            # now to write out each genotype for each gene
            genes = tsites[chrom][bp]

            for record in genes:
                gene  = record[0]
                etype = record[1]
                pval  = record[2]
                bh    = record[3]
                outline = [chrom, bp, etype, pval, bh, gene]
                genoExp = defaultdict(list) # genotype -> list[expr levels]
                for i in range(len(samples)):
                    sample   = samples[i]
                    genotype = genotypes[i]
                    ge_expr  = gexprs[sample][gene]
                    ge_expr  = round(ge_expr, 2)
                    entry    = f"{genotype}|{ge_expr}"
                    outline.append(entry)
                    if (ge_expr > 0):
                        genoExp[genotype].append(ge_expr)

                # get the avg normalized expression for each genotype
                avgs = []
                for genotype, exprs in genoExp.items():
                    avg = sum(exprs) / len(exprs)
                    avg = round(avg, 2)
                    avg = f"{len(exprs)}|{genotype}|{avg}"
                    avgs.append(avg)
                avgs = ", ".join(avgs)
                outline.append(avgs)

                # now write out the entry entry
                outline.append('\n')
                outline = '\t'.join(outline)
                outfh.write(outline)

    # close IOs
    fh.close()
    outfh.close()
    
    return 0

def main():
    """entry point to the program"""

    # get the inputs
    vcf, tgt, cnt = get_arguments()

    # get the targets
    get_targets(tgt)

    # get the expression levels for the target genes
    get_expression(cnt)

    # now construct the final result
    parse_vcf(vcf)


if __name__ == "__main__":
    main()