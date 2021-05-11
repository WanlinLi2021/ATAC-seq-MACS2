#!/bin/bash
#SBATCH --account=
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=30Gb
#SBATCH --time=5:00:00

module purge 2>/dev/null
module load mugqic/python/2.7.14
module load mugqic/MACS2/2.1.1.20160309

cd $SLURM_SUBMIT_DIR

macs2 callpeak \
--bdg \
--gsize mm \
--keep-dup all \
-f BEDPE \
--outdir ../../results/macs2/ \
-n 541_callpeak \
-t ../../results/filterLowQuality/541_sorted_rmMT_Dupli_lowQ.bam \
&> 541_callpeak.sh.log
