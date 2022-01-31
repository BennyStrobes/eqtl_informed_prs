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
		if len(data) != 4:
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

def get_num_eqtls(num_eqtls_file):
	f = open(num_eqtls_file)
	count = 0
	for line in f:
		line = line.rstrip()
		count = count + 1
	f.close()
	if count != 1:
		print('erroro')
		pdb.set_trace()
	return int(line)

def get_num_eqtls(num_eqtls_file, eqtl_counter):
	f = open(num_eqtls_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		data = line.split('\t')
		eqtl_counter[data[0]] = eqtl_counter[data[0]] + int(data[1])
	f.close()
	return eqtl_counter

def make_prs_beta_file_for_single_chromosome(chrom_num, gtex_tissues, gene_file, output_file, output_file2, output_file3, bivariate_cafeh_output_dir, version):
	num_tissues = len(gtex_tissues)
	# Initialize dicti to map from variant id to vector of length number of tissues (only have variant id if non-zero in one tissue)
	snp_counter = {}
	eqtl_counter = {}
	for tissue in gtex_tissues:
		snp_counter[tissue] = 0
		eqtl_counter[tissue] = 0
	dicti = {}
	dicti_weighted = {}

	data = np.loadtxt(gene_file, dtype=str, delimiter='\t')
	all_genes = data[1:, 0]
	all_chroms = data[1:, 1]
	indices = np.where(all_chroms==str(chrom_num))[0]
	genes = all_genes[indices]
	# Loop through genes
	for gene in genes:
		# get cafeh file for each gene
		cafeh_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name +'_' + version + '_' + gene + '_predicted_effects.txt'
		weighted_cafeh_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name +'_' + version + '_' + gene + '_weighted_predicted_effects.txt'
		raw_cafeh_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name +'_' + version + '_' + gene + '_pi.txt'
		num_eqtls_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name +'_' + version + '_' + gene + '_number_eqtl_components.txt'
		num_coloc_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name +'_' + version + '_' + gene + '_number_coloc_components.txt'
		# Extract number of eqtls identified by this tissue for this gene
		if os.path.exists(num_eqtls_file) == False:
			if os.path.exists(raw_cafeh_file) == True:
				print('assumption erorororooror')
				pdb.set_trace()
			continue
		# Keep track of number of components in each tissue
		eqtl_counter = get_num_eqtls(num_eqtls_file, eqtl_counter)
		snp_counter = get_num_eqtls(num_coloc_file, snp_counter)

		# check if file exists
		if os.path.exists(raw_cafeh_file) == False:
			continue
		# get effect sizes
		f = open(cafeh_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			effect_sizes = np.asarray(data[1:]).astype(float)
			if sum(effect_sizes != 0) == 0:
				continue
			snp_id = data[0]
			if snp_id not in dicti:
				dicti[snp_id] = np.zeros(num_tissues)
			for tissue_index in range(num_tissues):
				if dicti[snp_id][tissue_index] == 0.0:
					dicti[snp_id][tissue_index] = effect_sizes[tissue_index]
				else:
					dicti[snp_id][tissue_index] = np.mean([dicti[snp_id][tissue_index], effect_sizes[tissue_index]])
		f.close()
		# get effect sizes
		f = open(weighted_cafeh_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			effect_sizes = np.asarray(data[1:]).astype(float)
			if sum(effect_sizes != 0) == 0:
				continue
			snp_id = data[0]
			if snp_id not in dicti_weighted:
				dicti_weighted[snp_id] = np.zeros(num_tissues)
			for tissue_index in range(num_tissues):
				if dicti_weighted[snp_id][tissue_index] == 0.0:
					dicti_weighted[snp_id][tissue_index] = effect_sizes[tissue_index]
				else:
					dicti_weighted[snp_id][tissue_index] = np.mean([dicti_weighted[snp_id][tissue_index], effect_sizes[tissue_index]])
		f.close()


	t = open(output_file,'w')
	t.write('snp_name\t' + '\t'.join(gtex_tissues) + '\n')
	for snp_name in dicti.keys():
		t.write(snp_name + '\t' + '\t'.join(dicti[snp_name].astype(str)) + '\n')
	t.close()

	t = open(output_file2,'w')
	t.write('snp_name\t' + '\t'.join(gtex_tissues) + '\n')
	for snp_name in dicti_weighted.keys():
		t.write(snp_name + '\t' + '\t'.join(dicti_weighted[snp_name].astype(str)) + '\n')
	t.close()

	t = open(output_file3, 'w')
	t.write('tissue_name\tnum_components\tnum_eqtl_components\n')
	for tissue in gtex_tissues:
		t.write(tissue + '\t' + str(snp_counter[tissue]) + '\t' + str(eqtl_counter[tissue]) + '\n')
	t.close()



gtex_tissue_file = sys.argv[1]
trait_name = sys.argv[2]
processed_bivariate_cafeh_input_dir = sys.argv[3]
bivariate_cafeh_output_dir = sys.argv[4]
chrom_num = int(sys.argv[5])
version = sys.argv[6]

# Extract gtex tissues
gtex_tissues = extract_gtex_tissue_names(gtex_tissue_file)



gene_file = processed_bivariate_cafeh_input_dir + trait_name + '_processed_gene_list.txt'


# Make PRS beta file seperately for each chromosome
output_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_' + version + '_prs_beta_chrom_' + str(chrom_num) + '.txt'
output_file2 = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name +'_' + version + '_prs_weighted_beta_chrom_' + str(chrom_num) + '.txt'
output_file3 = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_' + version + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

make_prs_beta_file_for_single_chromosome(chrom_num, gtex_tissues, gene_file, output_file, output_file2, output_file3, bivariate_cafeh_output_dir, version)
