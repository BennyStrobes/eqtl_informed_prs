#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=50G                         # Memory total in MiB (for all cores)




cafeh_gene_list_file="$1"
processed_gtex_associations_dir="$2"
chrom_num="$3"
liftover_directory="$4"
eqtl_summary_stats_dir="$5"
gtex_tissue_file="$6"
genotype_reference_panel_dir="$7"

source ~/.bash_profile
python3 process_gtex_associations_for_cafeh_in_parallel.py $cafeh_gene_list_file $processed_gtex_associations_dir $chrom_num $eqtl_summary_stats_dir $gtex_tissue_file $genotype_reference_panel_dir


python3 liftover_cafeh_association_files_to_hg19.py $cafeh_gene_list_file $processed_gtex_associations_dir $liftover_directory $chrom_num
