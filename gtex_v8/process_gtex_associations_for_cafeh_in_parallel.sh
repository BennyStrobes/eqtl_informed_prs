#!/bin/bash -l

#SBATCH
#SBATCH --time=15:00:00
#SBATCH --partition=shared
#SBATCH --mem=5GB
#SBATCH --nodes=1



cafeh_gene_list="$1"
processed_gtex_associations_dir="$2"
job_number="$3"
num_jobs="$4"
liftover_directory="$5"



module load python/2.7-anaconda53
python process_gtex_associations_for_cafeh_in_parallel.py $cafeh_gene_list $processed_gtex_associations_dir $job_number $num_jobs

python liftover_cafeh_association_files_to_hg19.py $cafeh_gene_list $processed_gtex_associations_dir $liftover_directory $job_number $num_jobs
