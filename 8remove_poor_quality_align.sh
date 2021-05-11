#!/bin/bash
#SBATCH --account=
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=16Gb
#SBATCH --time=6:00:00

module purge 2>/dev/null
module load samtools/1.10
cd $SLURM_SUBMIT_DIR

samtools view \
-b \
-q 10 \
../../results/picard/541_sorted_rmMT_rmDupli.bam > ../../results/filterLowQuality/541_sorted_rmMT_Dupli_lowQ.bam

