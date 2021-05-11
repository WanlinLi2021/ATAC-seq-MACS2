#!/bin/bash
#SBATCH --account=
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=30Gb
#SBATCH --time=8:00:00

module purge 2>/dev/null
module load mugqic/java
module load StdEnv/2020
module load mugqic/picard/1.123
cd $SLURM_SUBMIT_DIR

java -jar /cvmfs/soft.mugqic/CentOS6/software/picard/picard-tools-1.123/CollectInsertSizeMetrics.jar \
I=../../results/picard/541_sorted_rmMT_rmDupli.bam \
O=../../results/taille_insert/541_insert_size_metrics.txt \
H=../../results/taille_insert/541_insert_size_histogram.pdf \
M=0.5 \
&> 541_insertSize.sh.log
