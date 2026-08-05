# README #

This repository is used to house the scripts used in the eQTL analysis pipeline for stickleback. The project is a joint effort from current and past members of the [Bell](https://alisonbelllab.web.illinois.edu/people/) and [Catchen](https://catchenlab.life.illinois.edu/) labs at the University of Illinois Urbana-Champaign.

The goal of the code found here is to identify nucleotide variants associated with transcriptomic differences between the common and white stickleback ecotypes.

First, the RNA-Seq data is processed through `STAR` and `HTSeq` to generate feature counts (see scripts in *bash_scripts*). The resulting counts are then normalized using `edgeR`.

On the genomic side of things, the genomic SNPs (from the *.vcf* file) were filtered to retain markers for individuals that had both RAD-Seq and RNA-Seq data. In addition, the distance the SNP was from the linkage marker determined whether the SNP will be retained. Out of the 166 samples with RAD data and the 70 with RNA-Seq data, 70 samples were retained. The filteration was performed with the `FilterVcfUsingLnkMap.py` script.

```
python3 FilterVcfUsingLnkMap.py --samples_list list_of_samples.txt --vcf populations.snps.vcf.gz --linkage_map linkage_map.csv --distance 100 --marker_coordiantes marker_coordinates.csv
```

An optional but recommended approach for further filtering is to use a minimum number of individuals required for a specific genotype. This was advice obtained from this [post](https://www.biostars.org/p/381662/), since `MatrixEQTL` can produce a tons of false positives due to rare genotypes. This filtering can be performed with `filter_vcf_genotype_counts.py`.
```
# filter the vcf file for variants where at least 5 individuals support each genotype
python3 filter_vcf_genotype_counts.py -v input.vcf.gz --count 5
```


The above analysis generates a vcf, which was converted into a matrix compatible with `MatrixEQTL` using another python script

```
python3 VcfToSNPsMatrix.py --vcf filtered_snps.vcf
```

Since not all the samples used in the `HTSeq` analysis have genomic data, the counts file was also filtered for the 70 remaining samples

```
python3 Subsample_htseq_counts.py --counts htseq_counts.txt --samples list_of_samples.txt
```

Lastly, `MatrixEQTL` requires the positions of the SNPs and Genes to distinguish *cis* from *trans* elements. These two files were created using either the *.gtf* file used in the RNA-Seq processing stage or vcf containing the genomic SNPs.

```
# get gene positions
python3 GTF2GenePos.py --gtf stickleback_v5_ensembl_genes.gtf.gz 

# get SNP positions
python3 VCF2SnpsPos.py -v filtered_snps.vcf
```

Altogether, `MatrixEQTL` uses 4 input files:

1. A SNPs matrix
2. Normalized featured counts
3. SNP positions
4. Gene positions

Once all the data is loaded, one can identify potentially associative *cis* and *trans* elements given a user defined distance.

Since there will be many eQTLs reported but not all may be related due to ecotypic differences (i.e., common vs white), there is an additional script that was written to filter for eQTLs for variants where the variants are either fixed between the two different ecotypes, or variants where they are either fixed between the two ecotypes or where one of the alleles is private in one of the ecotypes. The ecotypic information currently is hard coded in the python script but this will likely change when we obtain a large dataset of the two ecotypes.
```
python3 find_informative_sites.py --vcf ecotypes.vcf.gz --eqtls ematrix.results.tsv
```