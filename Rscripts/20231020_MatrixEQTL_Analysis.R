# set some general input variables
wd = "/Users/giovannimadrigal/Library/CloudStorage/Box-Box/USB_GIO/Research/Stickleback_eQTL/R_Analysis/20231020_Analysis/10Kb_Distance"
setwd(wd)

norm_counts_file = "gac_GE_Matrix.tsv"
snp_matrix_file = "Filtered_populations.snps_10000bp_Dist_SNPs_SNPsMatrix.tsv"

# setting some params for things I don't have
library(MatrixEQTL)

# following a guide from here: https://github.com/andreyshabalin/MatrixEQTL

# set the the inputs (just snps & gene expression - no covariates)
snps = SlicedData$new()
snps$fileDelimiter = '\t'
snps$fileOmitCharacters = 'NA' # missing vals
snps$fileSkipRows = 1 # skip sample names
snps$fileSkipColumns = 1 # skip snp IDs
snps$LoadFile(snp_matrix_file)

gene = SlicedData$new()
gene$fileDelimiter = '\t'
gene$fileOmitCharacters = "NA"
gene$fileSkipColumns = 1
gene$fileSkipRows = 1
gene$LoadFile(norm_counts_file)


# create the model using some testing params
output_file_name = paste(wd, "MatrixEQTL_Test.txt", sep = '/')
useModel = modelLINEAR # other options include modelANOVA, modelLINEAR_CROSS
pval_threshold = 1e-3
cisDist = 1e6
pval_threshold_cis = 1e-3

# now to include snp and gene positions
snpspos = "Filtered_populations.snps_10000bp_Dist_SNPs_snpspos.tsv"
genepos = "stickleback_v5_ensembl_genes_genepos_full.tsv"

snpspos = read.table(snpspos, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
genepos = read.table(genepos, header = TRUE, stringsAsFactors = FALSE, sep = '\t')


output_cis = "eQTL_Testing_cis"
me = Matrix_eQTL_main(snps = snps, gene = gene, cvrt = SlicedData$new(),
                        output_file_name = output_file_name, pvOutputThreshold = pval_threshold,
                        useModel = useModel, errorCovariance = numeric(), output_file_name.cis = output_cis,
                        verbose = TRUE, pvalue.hist = "qqplot", min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE, snpspos = snpspos, genepos = genepos,
                        pvOutputThreshold.cis = pval_threshold_cis, cisDist = cisDist)

View(me$cis$eqtls)












































