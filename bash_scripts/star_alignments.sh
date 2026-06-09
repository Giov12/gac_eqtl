#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --nodes=1 #nodes
#SBATCH --ntasks=24 #threads
#SBATCH --job-name=gac.star
#SBATCH --partition=catchenlab

# params
working_dir=/projects/catchenlab/gm33/gac_eQTL
reads_dir=${working_dir}/raw
outdir=${working_dir}/alignments
ref_dir=${working_dir}/assembly
gtf=${working_dir}/annotations/stickleback_v5_ensembl_genes.gtf
thr=${SLURM_NTASKS}


# start the alignment
cd $outdir

# going to align each sample individually
for sample in `ls -1 ${reads_dir}/*_R1_*`; do
    sample_name=$(basename $sample | cut -f 1-3 -d '_')
    fq1=${sample_name}_R1_001.fastq.gz
    fq2=${sample_name}_R2_001.fastq.gz
    STAR --runMode alignReads --genomeDir $ref_dir  \
        --sjdbGTFfile $gtf --sjdbOverhang 149  \
        --readFilesIn ${reads_dir}/${fq1} ${reads_dir}/${fq2} \
        --readFilesCommand zcat --runThreadN $thr  \
        --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $sample_name \
	    --outSAMunmapped Within --outBAMsortingThreadN $thr
done
