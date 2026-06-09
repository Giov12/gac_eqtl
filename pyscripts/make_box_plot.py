#!/bin/env python3

import gzip
import argparse
import os, sys
import random
import matplotlib.pyplot as plt
from   collections import defaultdict

chrom    = ''
bp       = ''
dataFile = ''
gene     = ''
type_    = ''
genoMap  = defaultdict(list)

def get_arguments() -> int:
    """function to get the target eQTL"""

    global chrom, bp, dataFile, gene

    parser = argparse.ArgumentParser(description="Draw a boxplot for the target eQTL")
    parser.add_argument("-c", "--chrom", required=True, type=str, help="chromosome name")
    parser.add_argument("-b", "--bp", required=True, type=str, help="basepair position")
    parser.add_argument("-g", "--gene", required=True, type=str, help="gene_id")
    parser.add_argument("-i", "--input", required=True, type=str, help="full expression, genotype file")

    args     = parser.parse_args()
    chrom    = args.chrom
    bp       = args.bp
    gene     = args.gene
    dataFile = args.input

    assert os.path.isfile(dataFile), f"Could not locate {dataFile}"
    assert bp.isdigit(), "--bp must be an integer"

    return 0

def load_data() -> int:
    """parse the tsv file for the target genes"""

    global dataFile, chrom, bp, genoMap, gene, type_

    fh = gzip.open(dataFile, "rt") if dataFile.endswith(".gz") else open(dataFile, 'r')
    fd = False

    for line in fh:
        if (len(line) == 0):
            continue
        fields = line.split('\t')
        if (fields[0] != chrom or fields[1] != bp or fields[5] != gene):
            continue
        fd    = True
        type_ = fields[2]
        for i in range(6, len(fields) - 2):
            entry = fields[i].split('|')
            geno  = entry[0]
            exp   = float(entry[1])
            genoMap[geno].append(exp)
        break

    fh.close()

    if (fd == False):
        msg = f"Could not find eQTL {chrom} {bp} for gene {gene} in {dataFile}"
        sys.exit(msg)

    return 0

def make_box_plot() -> int:
    """draw a box plot for the genotype expression map for the specific eqtl"""

    global chrom, bp, gene, genoMap

    genotypes   = sorted(genoMap.keys())
    expressions = [genoMap[geno] for geno in genotypes]

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.boxplot(
        expressions,
        positions=range(len(genotypes)),
        widths=0.4,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=2),
        boxprops=dict(facecolor="lightsteelblue", alpha=0.6),
        whiskerprops=dict(linewidth=1.2),
        capprops=dict(linewidth=1.2),
    )

    # Jittered points
    for i, values in enumerate(expressions):
        jitter = [random.uniform(-0.15, 0.15) + i for _ in values]
        ax.scatter(
            jitter, values,
            s=30, color="steelblue", alpha=0.7,
            edgecolors="white", linewidths=0.4,
            zorder=3,
        )

    ax.set_xticks(range(len(genotypes)))
    ax.set_xticklabels(genotypes, fontsize=12)
    ax.set_title(gene, fontsize=15)
    ax.set_xlabel( f"${type_}$ eQTL: {chrom} {bp}", fontsize=13)
    ax.set_ylabel("Normalized Expression", fontsize=13)
    ax.spines[["top", "right"]].set_visible(False)

    outfile = f"{gene}_{chrom}_{bp}_boxplot.png"
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)

    return 0

def main():
    """entry point to the pipeline"""

    # get the inputs
    get_arguments()

    # load the data
    load_data()

    # plot
    make_box_plot()

if __name__ == "__main__":
    main()