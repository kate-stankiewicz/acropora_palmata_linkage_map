#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N gemme_proc

cd $PBS_O_WORKDIR

conda deactivate

for phen in $(seq 1 12); do
	for spec in palm cerv; do
		Rscript manhattan_admixture_mapping.R ${phen} ${spec}
	done
done
