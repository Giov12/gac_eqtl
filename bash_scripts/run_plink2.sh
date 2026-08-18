#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1 #nodes
#SBATCH --ntasks=24 #threads
#SBATCH --mem=4G
#SBATCH --job-name=gac.plink2
#SBATCH --partition==<Partition-Name>

# params
base=/path/to/baseDir
working_dir=${base}/out_directory
vcf=${working_dir}/snp.vcf.gz
si=${working_dir}/sex_info.tsv
plink2=${base}/path/to/plink2
thr=${SLURM_NTASKS}

cd $working_dir

#
# sex_info.tsv is a 3 column .tsv file
#   
# column_1    column_2     Sex
# sample_name sample_name  1
# sample_name sample_name  1
# sample_name sample_name  1
#
# it contains all the same names in the vcf file
# for both column 1 and 2, and since all are male,
# they are encoded as 1
#

$plink2 --vcf $vcf --make-king square --pca --out gac.kinship.all \
        --threads $thr --update-sex $si --allow-extra-chr

