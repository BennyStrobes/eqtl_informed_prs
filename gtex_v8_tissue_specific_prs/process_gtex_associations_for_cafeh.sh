#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=60G                         # Memory total in MiB (for all cores)



processed_gtex_associations_dir="$1"
cafeh_gene_list_file="$2"
liftover_directory="$3"
eqtl_summary_stats_dir="$4"
gtex_tissue_file="$5"
genotype_reference_panel_dir="$6"

source ~/.bash_profile


for chrom_num in $(seq 1 22); do 
	sbatch process_gtex_associations_for_cafeh_in_parallel.sh $cafeh_gene_list_file $processed_gtex_associations_dir $chrom_num $liftover_directory $eqtl_summary_stats_dir $gtex_tissue_file $genotype_reference_panel_dir
done

