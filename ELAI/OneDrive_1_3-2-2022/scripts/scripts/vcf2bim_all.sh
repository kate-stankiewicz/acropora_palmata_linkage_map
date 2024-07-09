#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N v2b_a

cd $PBS_O_WORKDIR
path="/storage/home/khs18/scratch/genome_vcfs/chrom_files"

#also run same conversion on vcfs that are separated by chrom AND by species
for chr in $(cat chroms.txt); do
	for file in $(cat species.txt); do
		/gpfs/group/ibb3/default/tools/plink_v1.90/plink --vcf ${path}/${chr}/${file}_${chr}.v.vcf.recode.vcf --recode bimbam-1chr --aec --keep-allele-order --out ${path}/${chr}/${file}_${chr}_v_k
	done
done

