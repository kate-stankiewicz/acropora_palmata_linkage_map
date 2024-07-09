#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N c_id

cd $PBS_O_WORKDIR
conda activate base

#change sample name in vcf from affy ID to Baums ID
bcftools reheader -s newname.txt -o no_un_new_sampID_snpc_apalm_ref.vcf Apalm_ref_fixed_cc_snpchip.vcf
