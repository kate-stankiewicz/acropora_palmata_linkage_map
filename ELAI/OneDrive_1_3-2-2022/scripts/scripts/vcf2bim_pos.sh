#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N vcf2bim

cd $PBS_O_WORKDIR

#convert vcf to bimbam to get pos per chrom using vcfs per chrom (each including all 3 species)
for i in $(cat chrom_all_spec_vcfs.txt); do
/gpfs/group/ibb3/default/tools/plink_v1.90/plink --vcf $i --recode bimbam-1chr --aec --keep-allele-order --out ./all_pos/${i/.recode.vcf/}_k
done

 
