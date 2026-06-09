#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --nodes=1 #nodes
#SBATCH --ntasks=20 #threads
#SBATCH --job-name=gac.build.genome
#SBATCH --partition=catchenlab

# params
working_dir=/projects/catchenlab/gm33/gac_eQTL
outdir=${working_dir}/assembly
ref=${outdir}/stickleback_v5_assembly.fa
gtf=${working_dir}/annotations/stickleback_v5_ensembl_genes.gtf
thr=${SLURM_NTASKS}

# uncompress and will re-compress afterwards
gunzip ${ref}.gz

# start build
cd $outdir

# the gtf file created from gff3 using gffread
STAR --runThreadN $thr --runMode genomeGenerate --genomeDir ${outdir} \
    --genomeFastaFiles $ref --sjdbGTFfile $gtf --sjdbOverhang 149 \
    --genomeSAindexNbases 13 # based on previous run to avoid warning

# recompress
gzip $ref
