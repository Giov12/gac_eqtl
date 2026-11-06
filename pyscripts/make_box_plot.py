#!/bin/env python3

import gzip
import argparse
import os, sys
import random
import matplotlib.pyplot as plt
from   collections import defaultdict

rec_id   = ''
dataFile = ''
genoMap  = defaultdict(list)

def get_arguments() -> int:
    """function to get the target eQTL"""

    global rec_id, dataFile

    parser = argparse.ArgumentParser(description="Draw a boxplot for the target eQTL")
    parser.add_argument("-I", "--ID",    required=True, type=str, help="ID of record to plot")
    parser.add_argument("-i", "--input", required=True, type=str, help="full expression, genotype file")

    args     = parser.parse_args()
    rec_id   = args.ID
    dataFile = args.input

    assert os.path.isfile(dataFile), f"Could not locate {dataFile}"
    assert rec_id.isdigit(), "--ID is expected to be an integer"

    return 0

def load_data() ->tuple[str, str, str, str]:
    """parse the tsv file for the target genes"""

    global dataFile, rec_id, genoMap

    fh = gzip.open(dataFile, "rt") if dataFile.endswith(".gz") else open(dataFile, 'r')
    fd = False
    of = -1 # offset
    for line in fh:
        if (len(line) == 0):
            continue 
        if (of == -1):
            # set the offset
            if ("Allele.Info" in line):
                of = 2
            else:
                of = 1 
        fields = line.split('\t')
        if (fields[0] != rec_id):
            continue
        fd    = True
        chrom = fields[1]
        bp    = fields[2]
        type_ = fields[3]
        gene  = fields[6]
        for i in range(7, len(fields) - of):
            entry = fields[i].split('|')
            geno  = entry[0]
            exp   = float(entry[1])
            genoMap[geno].append(exp)
        break

    fh.close()

    if (fd == False):
        msg = f"Could not find eQTL {chrom} {bp} for gene {gene} in {dataFile}"
        sys.exit(msg)

    return (chrom, bp, gene, type_)

def make_box_plot() -> int:
    """draw a box plot for the genotype expression map for the specific eqtl"""

    global chrom, bp, gene, genoMap

    chrom, bp, gene, type_ = load_data()
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

    # plot
    make_box_plot()

if __name__ == "__main__":
    main()