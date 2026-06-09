# start off with setting up the directory and preparing input files
wd = "/Users/giovannimadrigal/Library/CloudStorage/Box-Box/USB_GIO/Research/Stickleback_eQTL/R_Analysis/20231019_Analysis"
setwd(wd)

# NOTE: txt file is really a tsv file
cnt_file = "/Users/giovannimadrigal/Library/CloudStorage/Box-Box/USB_GIO/Research/Stickleback_eQTL/htseq_counts/gac_rna_htseq_counts_with_header_subsampled.txt"

# libraries for the analysis
library(limma)
library(edgeR)

# okay, let's read in all the records and get the dimensions
gene_counts = read.delim(cnt_file, sep = '\t', 
                     row.names = 1, stringsAsFactors = FALSE, header = TRUE)
dim(gene_counts)
# [1] 15700    70

# I guess we can start creating the deg?
dge = DGEList(counts = gene_counts)

# it looks like a common practice is to filter generally lowly expressed genes
keep = filterByExpr(dge)
# probably need a group


# normalize
dge = calcNormFactors(dge)

# taking advice from this post to export / create normalized counts : https://www.biostars.org/p/9475236/
# TMM normalizes the library sizes to produce effective library sizes. 
# cpm values are counts normalized by the effective library sizes. 
# rpkm values are counts normalized by effective library sizes and by gene/feature length.

# no gene lengths in the obj, so using cpm
# NOTE: there is a mention to not use rpkm for these type of analyzes : https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
cpm_dge = cpm(y = dge)
cpm_dge

# now to write the output
dim(cpm_dge)
# [1] 15700    77

# write output
colnames(cpm_dge) = gsub(colnames(cpm_dge), pattern = "X", replacement = '')
write.table(cpm_dge, sep = '\t', row.names = TRUE, col.names = TRUE, 
            file = paste(wd, "gac_GE_Matrix.tsv", sep = '/'), quote = FALSE)



























