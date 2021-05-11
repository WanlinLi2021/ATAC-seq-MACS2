#!/bin/bash
#SBATCH --account= 
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=6Gb
#SBATCH --time=2:00:00

module purge 2>/dev/null
module load samtools/1.10
cd $SLURM_SUBMIT_DIR

samtools idxstats ../../results/bowtie2/541_sorted.bam | cut -f -1 | grep -v MT | xargs samtools view -b ../../results/bowtie2/541_sorted.bam > ../../results/removeMT/541_sorted_rmMT.bam
