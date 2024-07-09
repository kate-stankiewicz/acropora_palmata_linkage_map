#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N elai

cd $PBS_O_WORKDIR
path="/storage/home/khs18/scratch/genome_vcfs/chrom_files"


#run elai, loop through files for each chromosome

for chr in $(cat chroms.txt); do 
	for mg in 1 2 3; do
		/gpfs/group/ibb3/default/tools/elai/elai-lin -g ${path}/${chr}/palm_${chr}_v_k.recode.geno.txt -p 10 -g ${path}/${chr}/cerv_${chr}_v_k.recode.geno.txt -p 11 -g ${path}/${chr}/hybr_${chr}_v_k.recode.geno.txt -p 1 -pos ${path}/all_pos/all_${chr}_k.recode.pos.txt -o rr__${mg}_${chr}  -C 2 -c 10 -s 20 --exclude-nopos --exclude-maf 0.01 --exclude-miss1 -mg ${mg}
	done
done
