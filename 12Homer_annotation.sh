#!/bin/bash
#SBATCH --account=
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=3
#SBATCH --time=6:00:00

module purge 2>/dev/null
module load mugqic/homer/4.11

cd $SLURM_SUBMIT_DIR

annotatePeaks.pl \
../../results/macs3/541_callpeak_peaks.narrowPeak \
mm10 -gsize mm10 -cons -CpG \
-go ../../results/homer_annotatePeaks_Macs/541/gene_ontology \
-genomeOntology ../../results/homer_annotatePeaks_Macs/541/genome_ontology \
1> ../../results/homer_annotatePeaks_Macs/541_peaks_Macs.csv \
2> 541_macs_peaksAnnot.log
