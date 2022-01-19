#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1



processed_gtex_associations_dir="$1"
gene_info_file="$2"
liftover_directory="$3"


module load python/2.7-anaconda53

cafeh_gene_list=$processed_gtex_associations_dir"cafeh_gene_list.txt"
if false; then
python get_gtex_gene_list_for_cafeh.py $cafeh_gene_list $gene_info_file
fi


# Get summary stat association statistics for each gene
# Also liftover those sumstats from hg38 to hg19
num_jobs="20"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch process_gtex_associations_for_cafeh_in_parallel.sh $cafeh_gene_list $processed_gtex_associations_dir $job_number $num_jobs $liftover_directory
done
fi
