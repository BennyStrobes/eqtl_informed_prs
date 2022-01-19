import numpy as np 
import os
import sys
import pdb



def extract_ordered_genes_from_gene_list(cafeh_gene_list):
	f = open(cafeh_gene_list)
	genes = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes.append(data[0])
	return np.asarray(genes)

def association_file_to_bed_file(cafeh_hg38_file, temp_hg38_bed):
	f = open(cafeh_hg38_file)
	t = open(temp_hg38_bed,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 11:
			print('assumption eroror')
			pdb.set_trace()
		snp_id = data[3]
		snp_info = snp_id.split('_')
		chrom_num = snp_info[0]
		pos = snp_info[1]
		pos2 = int(pos) + 1
		t.write(chrom_num + '\t' + pos + '\t' + str(pos2) + '\n')
	f.close()
	t.close()

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg38ToHg19.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)


def recreate_cafeh_association_file_in_hg19_format(hg38_association_file, hg19_association_file, liftover_output_file, liftover_missing_file, hg38_to_hg19_mapping_file):
	missing = {}
	f = open(liftover_missing_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		missing[data[0] + '_' + data[1]] = 1
	f.close()
	f = open(hg38_association_file)
	t = open(hg19_association_file,'w')
	t_mapping = open(hg38_to_hg19_mapping_file,'w')
	g = open(liftover_output_file)
	t_mapping.write('hg38_snp_id\thg19_snp_id\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		snp_id = data[3]
		snp_info = snp_id.split('_')
		snp_key = snp_info[0] + '_' + snp_info[1]
		if snp_key in missing:
			t_mapping.write(snp_id + '\tNA\n')
			continue
		hg19_info = next(g).rstrip().split('\t')
		new_snp_id = hg19_info[0] + '_' + hg19_info[1] + '_' + snp_info[2] + '_' + snp_info[3]
		t.write(','.join(data[:3]) + ',' + new_snp_id + ',' + ','.join(data[4:]) + '\n')
		t_mapping.write(snp_id + '\t' + new_snp_id + '\n')
	f.close()
	g.close()
	t.close()
	t_mapping.close()

def get_num_lines(file_name):
	counter = 0
	f = open(file_name)
	for line in f:
		counter = counter + 1
	f.close()
	return counter

def error_checking(liftover_output_file, hg19_association_file):
	num_lines_lift = get_num_lines(liftover_output_file)
	num_lines_association = get_num_lines(hg19_association_file) - 1
	if num_lines_association != num_lines_lift:
		print('assumption erororo')



# For parallelization purposes
def parallelization_start_and_end(num_tasks, job_number, total_jobs):
    tasks_per_job = (num_tasks/total_jobs) + 1
    start_task = job_number*tasks_per_job
    end_task = (job_number + 1)*tasks_per_job -1 
    return start_task, end_task






cafeh_gene_list = sys.argv[1]
processed_gtex_associations_dir = sys.argv[2]
liftover_directory = sys.argv[3]
job_number = int(sys.argv[4])
total_jobs = int(sys.argv[5])



ordered_genes = extract_ordered_genes_from_gene_list(cafeh_gene_list)

#For parallelization purposes
start_number, end_number = parallelization_start_and_end(len(ordered_genes), job_number, total_jobs)

# Loop through genes
for gene_index, gene in enumerate(ordered_genes):
	# Skip genes not in this parallelization run
	if gene_index < start_number or gene_index > end_number:
		continue
	# Cafeh file in hg38 format
	cafeh_hg38_file = processed_gtex_associations_dir + gene + '_associations.csv'
	# First convert hg38  association to temporary bed format
	temp_hg38_bed = processed_gtex_associations_dir + gene + '_temp_bed_hg38.txt'
	association_file_to_bed_file(cafeh_hg38_file, temp_hg38_bed)

	# Run liftover
	liftover_output_file = processed_gtex_associations_dir + gene + '_liftover_output.txt'
	liftover_missing_file = processed_gtex_associations_dir + gene + '_liftover_missing.txt'
	run_liftover(temp_hg38_bed, liftover_output_file, liftover_missing_file, liftover_directory)


	# Cafeh association file in hg19 format
	cafeh_hg19_file = processed_gtex_associations_dir + gene + '_hg19_associations.csv'
	hg38_to_hg19_mapping_file = processed_gtex_associations_dir + gene + '_hg38_to_hg19_snp_mapping.txt'
	recreate_cafeh_association_file_in_hg19_format(cafeh_hg38_file, cafeh_hg19_file, liftover_output_file, liftover_missing_file, hg38_to_hg19_mapping_file)

	# Quick error checking
	error_checking(liftover_output_file, cafeh_hg19_file)


	# Remove temporary files
	os.system('rm ' + liftover_output_file)
	os.system('rm ' + liftover_missing_file)
	os.system('rm ' + temp_hg38_bed)

