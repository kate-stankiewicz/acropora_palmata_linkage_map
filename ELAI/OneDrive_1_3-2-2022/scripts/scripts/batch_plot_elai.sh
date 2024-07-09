#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=35gb
#PBS -j oe
#PBS -N plots

cd $PBS_O_WORKDIR

conda deactivate
module purge
module load r/3.4

for chrom in $(cat chroms.txt); do
	for mg in mg1 mg2 mg5; do
		Rscript elai_plots.R $chrom $mg
	done
done
