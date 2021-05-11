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

java -jar /cvmfs/soft.mugqic/CentOS6/software/picard/picard-tools-1.123/MarkDuplicates.jar \
INPUT=../../results/removeMT/541_sorted_rmMT.bam \
OUTPUT=../../results/picard/541_sorted_rmMT_rmDupli.bam \
METRICS_FILE=../../results/picard/541_picard_metrics.txt \
REMOVE_DUPLICATES=true \
&> 541_picard.sh.log
