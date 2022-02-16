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

def get_number_of_eqtls_from_p_active(p_active_file, version, eqtl_counter):
	p_active_raw = np.loadtxt(p_active_file, dtype=str,delimiter='\t')
	p_active_raw = p_active_raw[:-1,:]
	tissue_names = p_active_raw[:,0]
	p_active = p_active_raw[:,1:].astype(float)
	version_info = version.split('_')
	p_active_thresh = float(version_info[-1])
	method = '_'.join(version_info[:-1])
	if method == 'all_tissues':
		counts = np.sum(p_active > p_active_thresh,axis=1)
		if len(counts) != len(tissue_names):
			print('assumtpion eroror')
			pdb.set_trace()
		for itera, count in enumerate(counts):
			tissue = tissue_names[itera]
			eqtl_counter[tissue] = eqtl_counter[tissue] + count
	elif method == 'top_tissue':
		num_comp = p_active.shape[1]
		for comp_num in range(num_comp):
			top_tissue_indices = np.where(p_active[:,comp_num] == np.max(p_active[:,comp_num]))[0]
			for tissue_index in top_tissue_indices:
				if p_active[tissue_index, comp_num] > p_active_thresh:
					tissue = tissue_names[tissue_index]
					eqtl_counter[tissue] = eqtl_counter[tissue] + 1
	else:
		print('assumption eroror')
		pdb.set_trace()

	return eqtl_counter

def get_number_of_eqtls_from_coloc_pp_file(coloc_pp_file, thresh, eqtl_counter, coloc_counter):
	f = open(coloc_pp_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue = data[0]
		pps = np.asarray(data[1:]).astype(float)
		pph1 = pps[1]
		pph3 = pps[3]
		pph4 = pps[4]
		if pph4 > thresh:
			coloc_counter[tissue] = coloc_counter[tissue] + 1
		if (pph4 + pph3 + pph1) > thresh:
			eqtl_counter[tissue] = eqtl_counter[tissue] + 1
	f.close()
	return eqtl_counter, coloc_counter


def make_prs_beta_file_for_single_chromosome(trait_name, chrom_num, gtex_tissues, gene_file, output_file, output_file2, bivariate_cafeh_output_dir, coloc_threshy, version):
	num_tissues = len(gtex_tissues)
	# Initialize dicti to map from variant id to vector of length number of tissues (only have variant id if non-zero in one tissue)
	eqtl_counter = {}
	coloc_counter = {}
	for tissue in gtex_tissues:
		eqtl_counter[tissue] = 0
		coloc_counter[tissue] = 0
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
		if version == 'adaptive':
			predicted_effects_file = bivariate_cafeh_output_dir + trait_name + '_' + gene + '_adaptive_coloc_' + coloc_threshy + '_predicted_effect_sizes.txt'
			coloc_pp_file = bivariate_cafeh_output_dir + trait_name + '_' + gene + '_adaptive_coloc_' + 'posterior_probabilities.txt'
		elif version == 'standard':
			predicted_effects_file = bivariate_cafeh_output_dir + trait_name + '_' + gene + '_coloc_' + coloc_threshy + '_predicted_effect_sizes.txt'
			coloc_pp_file = bivariate_cafeh_output_dir + trait_name + '_' + gene + '_coloc_' + 'posterior_probabilities.txt'			


		eqtl_counter, coloc_counter = get_number_of_eqtls_from_coloc_pp_file(coloc_pp_file, float(coloc_threshy), eqtl_counter, coloc_counter)

		# Extract number of eqtls identified by this tissue for this gene
		if os.path.exists(predicted_effects_file) == False:
			continue

		# get effect sizes
		f = open(predicted_effects_file)
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


	t = open(output_file,'w')
	t.write('snp_name\t' + '\t'.join(gtex_tissues) + '\n')
	for snp_name in dicti.keys():
		t.write(snp_name + '\t' + '\t'.join(dicti[snp_name].astype(str)) + '\n')
	t.close()

	t = open(output_file2, 'w')
	t.write('tissue_name\tnum_coloc\tnum_eqtls\n')
	for tissue in gtex_tissues:
		t.write(tissue + '\t' + str(coloc_counter[tissue]) + '\t' + str(eqtl_counter[tissue]) + '\n')
	t.close()



gtex_tissue_file = sys.argv[1]
trait_name = sys.argv[2]
processed_bivariate_cafeh_input_dir = sys.argv[3]
bivariate_cafeh_output_dir = sys.argv[4]
chrom_num = int(sys.argv[5])

# Extract gtex tissues
gtex_tissues = extract_gtex_tissue_names(gtex_tissue_file)



gene_file = processed_bivariate_cafeh_input_dir + trait_name + '_processed_gene_list.txt'



coloc_thresholds = [.5, .7, .9, .95, .99]


for coloc_threshold in coloc_thresholds:

	# Make PRS beta file seperately for each chromosome
	output_file = bivariate_cafeh_output_dir + 'adaptive_prior_coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_prs_beta_chrom_' + str(chrom_num) + '.txt'
	output_file2 = bivariate_cafeh_output_dir + 'adaptive_prior_coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'
	make_prs_beta_file_for_single_chromosome(trait_name,chrom_num, gtex_tissues, gene_file, output_file, output_file2, bivariate_cafeh_output_dir,str(coloc_threshold), 'adaptive')


	# Make PRS beta file seperately for each chromosome
	output_file = bivariate_cafeh_output_dir + 'coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_prs_beta_chrom_' + str(chrom_num) + '.txt'
	output_file2 = bivariate_cafeh_output_dir + 'coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'
	make_prs_beta_file_for_single_chromosome(trait_name,chrom_num, gtex_tissues, gene_file, output_file, output_file2, bivariate_cafeh_output_dir,str(coloc_threshold), 'standard')

