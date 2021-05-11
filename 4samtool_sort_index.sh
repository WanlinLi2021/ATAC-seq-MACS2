#!/bin/bash
#SBATCH --account=
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=16Gb
#SBATCH --time=6:00:00

module purge 2>/dev/null
module load samtools/1.10
cd $SLURM_SUBMIT_DIR

samtools sort \
-o ../../results/bowtie2/541_sorted.bam \
../../results/bowtie2/541.bam \
&> samtools_sort_541.sh.log

samtools index ../../results/bowtie2/541_sorted.bam \
2>> samtools_sort_541.sh.log

