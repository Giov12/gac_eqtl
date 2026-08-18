#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --nodes=1 #nodes
#SBATCH --ntasks=48 #threads
#SBATCH --mem=12G
#SBATCH --job-name=gac.subread
#SBATCH --partition=<Partition-Name>

# params
featurecounts=/path/to/featureCounts
working_dir=/path/to/parentDir
outdir=${working_dir}/gene.counts
samples_file=${outdir}/samples.txt # 1 column list of alignment files (*.bam)
gtf=${working_dir}/annotations/assemblyV113/Gasterosteus_aculeatus.UGA_version5.113.fixMT.gtf
thr=${SLURM_NTASKS}

# start the alignment
cd $outdir

samples=$(cat $samples_file | tr '\n' ' ')

$featurecounts -p --countReadPairs -C -o gac.counts.txt -t gene -T $thr -a $gtf $samples

