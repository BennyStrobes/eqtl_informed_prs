import numpy as np 
import os
import sys
import pdb







def extract_gtex_tissue_names(gtex_tissue_file):
	f = open(gtex_tissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 2:
			print('assumption eroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()

	return np.asarray(arr)

def extract_num_components_from_raw_cafeh_file(raw_cafeh_file):
	f = open(raw_cafeh_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		dicti[data[3]] = 1
	f.close()
	return len(dicti)



def make_prs_beta_file_for_single_chromosome(chrom_num, gtex_tissues, gene_files, output_file, output_file2):
	num_tissues = len(gtex_tissues)
	# Initialize dicti to map from variant id to vector of length number of tissues (only have variant id if non-zero in one tissue)
	snp_counter = {}
	for tissue in gtex_tissues:
		snp_counter[tissue] = 0
	dicti = {}
	# loop through tissues
	for tissue_index, gene_file in enumerate(gene_files):
		tissue_name = gtex_tissues[tissue_index]
		print(tissue_name)
		data = np.loadtxt(gene_file, dtype=str, delimiter='\t')
		all_genes = data[1:, 0]
		all_chroms = data[1:, 1]
		indices = np.where(all_chroms==str(chrom_num))[0]
		genes = all_genes[indices]
		# Loop through genes
		for gene in genes:
			# get cafeh file for each gene
			cafeh_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_' + tissue_name + '_' + gene + '_predicted_effects.txt'
			raw_cafeh_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_' + tissue_name + '_' + gene + '_coloc_snps_fixed_pi_summary_table.txt'
			# check if file exists
			if os.path.exists(raw_cafeh_file) == False:
				continue
			# Keep track of number of components in each tissue
			num_snps = extract_num_components_from_raw_cafeh_file(raw_cafeh_file)
			snp_counter[tissue_name] = snp_counter[tissue_name] + num_snps
			# get effect sizes
			f = open(cafeh_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				if data[1] == '0.0':
					continue
				snp_id = data[0]
				effect_size = float(data[1])
				if snp_id not in dicti:
					dicti[snp_id] = np.zeros(num_tissues)
				if dicti[snp_id][tissue_index] == 0.0:
					dicti[snp_id][tissue_index] = effect_size
				else:
					dicti[snp_id][tissue_index] = np.mean([dicti[snp_id][tissue_index], effect_size])
			f.close()
	t = open(output_file,'w')
	t.write('snp_name\t' + '\t'.join(gtex_tissues) + '\n')
	for snp_name in dicti.keys():
		t.write(snp_name + '\t' + '\t'.join(dicti[snp_name].astype(str)) + '\n')
	t.close()
	t = open(output_file2, 'w')
	t.write('tissue_name\tnum_components\n')
	for tissue in gtex_tissues:
		t.write(tissue + '\t' + str(snp_counter[tissue]) + '\n')
	t.close()



gtex_tissue_file = sys.argv[1]
trait_name = sys.argv[2]
processed_bivariate_cafeh_input_dir = sys.argv[3]
bivariate_cafeh_output_dir = sys.argv[4]
chrom_num = int(sys.argv[5])

# Extract gtex tissues
gtex_tissues = extract_gtex_tissue_names(gtex_tissue_file)

gene_files = []
for gtex_tissue in gtex_tissues:
	gene_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gtex_tissue + '_processed_gene_list.txt'
	gene_files.append(gene_file)
gene_files = np.asarray(gene_files)

# Make PRS beta file seperately for each chromosome
print(chrom_num)
output_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_prs_beta_chrom_' + str(chrom_num) + '.txt'
output_file2 = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

make_prs_beta_file_for_single_chromosome(chrom_num, gtex_tissues, gene_files, output_file, output_file2)


