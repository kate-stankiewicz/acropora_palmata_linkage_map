#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N sub_spec

cd $PBS_O_WORKDIR
conda activate base

#set up directories
for dir in $(cat chroms.txt); do
mkdir -p chrom_files/${dir}
done

#try vcftools instead of bcftools to get separate file for each species per chr
for i in $(cat chroms.txt); do
       for s in $(cat species.txt); do
                /gpfs/group/ibb3/default/tools/vcftools --vcf ./chrom_files/${i}.recode.vcf --keep ${s}_IDs.txt --recode --out ./chrom_files/${i}/${s}_${i}.v.vcf
        done
done
