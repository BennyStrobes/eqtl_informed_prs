#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)




pseudotissue_name="$1"
composit_tissue_string="$2"
pseudotissue_eqtl_dir="$3"


echo $pseudotissue_name"\t"$composit_tissue_string

for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	python3 run_eqtl_meta_analysis_for_single_pseudotissue_chromosome.py $pseudotissue_name $composit_tissue_string $pseudotissue_eqtl_dir $chrom_num
done
