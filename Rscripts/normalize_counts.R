# start off with setting up the directory and preparing input files
wd = "path_to_results_dir"
setwd(wd)

#
# read in a tsv file
# NOTE: This file is a result from the software featureCounts
# this file also assumes the first row has columns with the sample names
# instead of the alignment files used
#

cnt_file = "features_counts.tsv"

# library for the analysis
library(edgeR)

# okay, let's read in all the records and get the dimensions
gene_counts = read.delim(cnt_file, sep = '\t', 
                         row.names = 1, stringsAsFactors = FALSE, header = TRUE)
dim(gene_counts)
# [1] GeneCount    SampleCount

# Create the deg list
dge = DGEList(counts = gene_counts)

# Calculate % reads per gene
percent_reads <- (dge$counts / colSums(dge$counts)) * 100

# examining percentage distribution of aligned reads (these will actually be smaller when considering all reads)
percent_reads[apply(percent_reads, 1, function(r) any (r > 2.0)),]

# following the same filtering as previous work from the Bell lab
keep <- rowSums(cpm(dge) > 0.5) >= 4  # Keep genes with >1 CPM in at least 2 samples
dge  <- dge[keep, , keep.lib.sizes = FALSE]


# TMM normalizes the library sizes to produce effective library sizes. 
# normalize
dge = calcNormFactors(dge, method = "TMM")

norm_cpm <- cpm(dge, log = FALSE)  # Raw CPM
log2_cpm <- cpm(dge, log = TRUE)   # Log2-transformed CPM (better for visualization)

# now to write the output
dim(norm_cpm)
# [1] FeatureCount    SampleCount

dim(log2_cpm)
# [1] FeatureCount    SampleCount

# write output
colnames(norm_cpm) = gsub(colnames(norm_cpm), pattern = "X", replacement = '')
write.table(norm_cpm, sep = '\t', row.names = TRUE, col.names = TRUE, 
            file = paste(wd, "gac.norm.counts.raw.tsv", sep = '/'), quote = FALSE)

colnames(log2_cpm) = gsub(colnames(log2_cpm), pattern = "X", replacement = '')
write.table(log2_cpm, sep = '\t', row.names = TRUE, col.names = TRUE, 
            file = paste(wd, "gac.norm.counts.log.tsv", sep = '/'), quote = FALSE)





















