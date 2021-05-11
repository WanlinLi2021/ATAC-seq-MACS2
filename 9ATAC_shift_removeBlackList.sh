#!/bin/bash
#SBATCH --account=
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=30Gb
#SBATCH --time=8:00:00

module purge 2>/dev/null
module load mugqic/deepTools/3.5.0
cd $SLURM_SUBMIT_DIR

alignmentSieve -b ../../results/filterLowQuality/541_sorted_rmMT_Dupli_lowQ.bam \
--BED \
--ATACshift \
--blackListFileName mm10.blacklist.bed.gz \
-o ../../results/shiftAndRmblackList/541_bienFilter.bedpe \
--filterMetrics ../../results/shiftAndRmblackList/541_mstrics_log.txt
