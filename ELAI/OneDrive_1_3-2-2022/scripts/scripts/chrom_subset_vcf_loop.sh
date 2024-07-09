#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N sub_chromL

cd $PBS_O_WORKDIR

#step 1: subset the vcf by chrom
dir="chrom_files"
mkdir -p $dir
 
for i in $(cat chroms.txt); do
	/gpfs/group/ibb3/default/tools/vcftools --vcf ./clonecorrect_biallelic_snp.vcf --snps ${i}_snpIDS.txt --recode --recode-INFO-all --out ./${dir}/${i}
done
