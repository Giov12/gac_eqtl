#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --nodes=1 #nodes
#SBATCH --ntasks=4 #threads
#SBATCH --job-name=gac.htseq-counts
#SBATCH --partition=catchenlab

# params
working_dir=/projects/catchenlab/gm33/gac_eQTL
gtf=${working_dir}/annotations/stickleback_v5_ensembl_genes.gtf
bams=$(ls -1 ${working_dir}/alignments/*.bam | tr '\n' ' ' | sed 's/ $//g')
outdir=${working_dir}/htseq_counts


# python interpreter & setting path
unset PYTHONPATH

module load python/2
PYTHONPATH=/usr/local/python/2.7.15/bin/python2
export PYTHONPATH

# run the software
# format = bam
# stranded = no
htseq-count -f bam -s no $bams $gtf > ${outdir}/gac_rna_htseq_counts.txt
