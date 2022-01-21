import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import rnaseqnorm
import pandas as pd
from sklearn.decomposition import PCA


# We are going to limit analysis to genes tested in all tissues
def get_genes_tested_in_all_tissues(tissues, gtex_expression_dir):
	# Initialize dictionaries to keep track of which genes were in all tissues
	gene_counts = {}
	valid_genes = {}
	# For each tissue, keep track fo which genes were used
	for tissue in tissues:
		gene_file_name = gtex_expression_dir + tissue + '.tpm.txt'
		head_count = 0
		f = open(gene_file_name)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[0]
			#tpms = np.asarray(data[1:]).astype(float)
			#if sum(tpms > .1)/len(tpms) < .2:
			#	continue
			# Add gene to dictionary if not observed
			if ensamble_id not in gene_counts:
				gene_counts[ensamble_id] = 0
			# Keep track of how many tissues this gene was observed in 
			gene_counts[ensamble_id] = gene_counts[ensamble_id] + 1
		f.close()
	# Loop through all observed ensamble ids
	# Only take genes that are expressed in all tissues
	for ensamble_id in gene_counts.keys():
		if gene_counts[ensamble_id] == len(tissues):
			valid_genes[ensamble_id] = 1
	return valid_genes

def get_pseudotissue_sample_names(file_name):
	arr = []
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	f.close()
	return np.asarray(arr)


def filter_lowly_expressed_genes(tpm_matrix, ordered_genes):
	valid_indices = []
	for gene_number in range(len(ordered_genes)):
		if sum(tpm_matrix[:,gene_number] >= .1)/len(tpm_matrix[:,gene_number]) >= .2:
			valid_indices.append(gene_number)
	valid_indices = np.asarray(valid_indices)
	return tpm_matrix[:,valid_indices], np.asarray(ordered_genes)[valid_indices]

# Generate TPM expression matrix
def generate_tpm_expression_matrix(tissues, ordered_sample_names, genes_tested_in_all_tissues, gtex_tpm_dir, output_file):
	# Generate mapping from sample_name to row number
	sample_name_to_row_number = {}
	for index, sample_name in enumerate(ordered_sample_names):
		sample_name_to_row_number[sample_name] = index
	# Generate mapping from gene to column number
	ordered_genes = sorted(genes_tested_in_all_tissues.keys())
	gene_to_column_number = {}
	for index, gene in enumerate(ordered_genes):
		gene_to_column_number[gene] = index
	# Num samples and number of genes
	num_samples = len(sample_name_to_row_number)
	num_genes = len(ordered_genes)
	# Initialize tpm matrix
	tpm_matrix = np.zeros((num_samples, num_genes))
	counter = 0
	used = {}
	# Loop through each tissue and fill in the tpm matrix
	for tissue_index, tissue in enumerate(tissues):
		# Stream tpm file in this tissue
		tpm_file = gtex_tpm_dir + tissue + '.tpm.txt'
		head_count = 0
		f = open(tpm_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Header
			if head_count == 0:
				head_count = head_count + 1
				# Get ordered list of sample names
				sample_names = []
				for sample_name_temp in data[1:]:
					sample_name = tissue + ':' + sample_name_temp
					sample_names.append(sample_name)
				continue
			# Get relevent fields from line
			ensamble_id = data[0]
			# Skip genes we are not interested in
			if ensamble_id not in gene_to_column_number:
				continue
			tpm_counts = np.asarray(data[1:]).astype(float)
			# Loop through samples and add sample/tpm count to tpm_matrix
			for index, sample_name in enumerate(sample_names):
				# Ignore samples that aren't in our list
				if sample_name not in sample_name_to_row_number:
					continue
				# Get row corresponding to this smample
				row_index = sample_name_to_row_number[sample_name]
				# Get column corresponding to this gene
				column_index = gene_to_column_number[ensamble_id]

				#ele_name = str(row_index) + '_' + str(column_index)

				# Add tpm count to matrix
				tpm_matrix[row_index, column_index] = tpm_counts[index]
				counter = counter + 1
		f.close()
	# Remove genes with zero expression across all samples
	filtered_tpm_matrix, filtered_orderd_genes = filter_lowly_expressed_genes(tpm_matrix, ordered_genes)
	print(tpm_matrix.shape)
	print(filtered_tpm_matrix.shape)
	#log_filtered_tpm_matrix = np.log2(filtered_tpm_matrix + 1.0)
	# Print to output file
	t = open(output_file, 'w')
	# print header
	t.write('GeneId\t' + '\t'.join(filtered_orderd_genes) + '\n')
	for sample_num, sample_name in enumerate(ordered_sample_names):
		tpm_expr = filtered_tpm_matrix[sample_num,:].astype(str)
		t.write(sample_name + '\t' + '\t'.join(tpm_expr) + '\n')
	t.close()

def standardize_expression(tpm_expression_matrix_file, standardized_tpm_expression_matrix_file):
	tpm_full = np.loadtxt(tpm_expression_matrix_file, dtype=str,delimiter='\t')
	tpm = tpm_full[1:,1:].astype(float)
	samples = tpm_full[1:,0]
	genes = tpm_full[0,1:]
	# Quantile normalize the samples
	df = pd.DataFrame(np.transpose(tpm))

	temp_out = rnaseqnorm.normalize_quantiles(df)
	norm_df = rnaseqnorm.inverse_normal_transform(temp_out)
	standardized_tpm = np.transpose(np.asarray(norm_df))

	# Print to output file
	t = open(standardized_tpm_expression_matrix_file, 'w')
	# print header
	t.write('SampleId\t' + '\t'.join(samples) + '\n')
	for gene_num, gene_name in enumerate(genes):
		#expr = tpm_quantile_normalized[sample_num, :].astype(str)
		###
		expr = standardized_tpm[:, gene_num].astype(str)
		###
		t.write(gene_name + '\t' + '\t'.join(expr) + '\n')
	t.close()

def generate_expression_pcs(expression_file, pc_file, num_expression_pcs):
	# Load in expression data
	expr_full = np.transpose(np.loadtxt(expression_file, dtype=str,delimiter='\t'))
	expr = expr_full[1:,1:].astype(float)
	samples = expr_full[1:,0]
	genes = expr_full[0,1:]

	# Compute gene expression pcs on expression data
	#uuu, sss, vh = np.linalg.svd(np.transpose(expr))
	#expr_pc_loadings = np.transpose(vh)[:,:num_expression_pcs]

	# This way of computing pcs is faster
	_pca = PCA(n_components=num_expression_pcs, svd_solver='arpack')
	expr_pc_loadings = _pca.fit_transform(expr)

	# print to output file 
	pc_names = []
	for pc_num in range(num_expression_pcs):
		pc_names.append('PC' + str(pc_num))
	t = open(pc_file, 'w')
	t.write('Sample_name\t' + '\t'.join(pc_names) + '\n')
	for sample_num, sample_name in enumerate(samples):
		t.write(sample_name + '\t' + '\t'.join(expr_pc_loadings[sample_num,:].astype(str)) + '\n')
	t.close()

def create_mapping_from_gene_name_to_gene_info(gene_annotation_file):
	f = open(gene_annotation_file)
	mapping = {}
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		if len(data) != 9:
			print('assumption eroror')
			pdb.set_trace()
		if data[2] != 'gene':
			continue
		ensamble_id = 'null'
		gene_type = 'null'
		gene_info = data[8].split(';')
		for info in gene_info:
			if info.startswith('gene_id'):
				ensamble_id = info.split('"')[1]
			elif info.startswith(' gene_type'):
				gene_type = info.split('"')[1]
		if ensamble_id == 'null' or gene_type == 'null':
			print('assumption eroror')
			pdb.set_trace()
		gene_chrom_num = data[0]
		gene_strand = data[6]
		if gene_strand == '+':
			tss = data[3]
		elif gene_strand == '-':
			tss = data[4]
		else:
			print('assumption error')


		# Add to info
		if ensamble_id not in mapping:
			mapping[ensamble_id] = (gene_type, gene_chrom_num, gene_strand, tss)
		else:
			if mapping[ensamble_id][0] != gene_type:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][1] != gene_chrom_num:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][2] != gene_strand:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][3] != tss:
				print('assumption eroror')
				pdb.set_trace()
	f.close()
	return mapping


def filter_expression_matrix(standardized_tpm_expression_matrix_file, standardized_tpm_expression_matrix_filtered_genes_file, gene_name_to_gene_info):
	valid_chromosomes = {}
	for chrom_num in range(1,23):
		valid_chromosomes['chr' + str(chrom_num)] = 1
	valid_gene_types = {}
	valid_gene_types['protein_coding'] = 1
	valid_gene_types['lincRNA'] = 1

	f = open(standardized_tpm_expression_matrix_file)
	t = open(standardized_tpm_expression_matrix_filtered_genes_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + 'gene_chrom_num\tgene_tss\tgene_strand\t' + '\t'.join(data[1:]) + '\n')
			continue
		ensamble_id = data[0]
		gene_info = gene_name_to_gene_info[ensamble_id]
		gene_type = gene_info[0]
		gene_chrom_num = gene_info[1]
		gene_strand = gene_info[2]
		gene_tss = gene_info[3]

		if gene_type not in valid_gene_types:
			continue
		if gene_chrom_num not in valid_chromosomes:
			continue
		# Passed filtering
		t.write(data[0] + '\t' + gene_chrom_num + '\t' + gene_tss + '\t' + gene_strand + '\t' + '\t'.join(data[1:]) + '\n')
	f.close()
	t.close()


pseudotissue = sys.argv[1]
composit_tissue_string = sys.argv[2]
pseudotissue_sample_names_dir = sys.argv[3]
gtex_tpm_expression_dir = sys.argv[4]
gene_annotation_file = sys.argv[5]
pseudotissue_expression_dir = sys.argv[6]




# Create mapping from gene name to gene class
gene_name_to_gene_info = create_mapping_from_gene_name_to_gene_info(gene_annotation_file)

# Get array of composit tissues
composit_tissues = np.asarray(composit_tissue_string.split(','))

# Get pseudotissue sample names
pseudotissue_ordered_sample_names = get_pseudotissue_sample_names(pseudotissue_sample_names_dir + pseudotissue + '_sample_names.txt')

# We are going to limit analysis to genes tested in all tissues
genes_tested_in_all_tissues = get_genes_tested_in_all_tissues(composit_tissues, gtex_tpm_expression_dir)



# Generate TPM expression matrix
tpm_expression_matrix_file = pseudotissue_expression_dir + pseudotissue + '_tpm.txt'
generate_tpm_expression_matrix(composit_tissues, pseudotissue_ordered_sample_names, genes_tested_in_all_tissues, gtex_tpm_expression_dir, tpm_expression_matrix_file)


# Quantile normalize and standardize TPM expression matrix
standardized_tpm_expression_matrix_file = pseudotissue_expression_dir + pseudotissue + '_normalized_expression.txt'
standardize_expression(tpm_expression_matrix_file, standardized_tpm_expression_matrix_file)


# Filter genes based on autosomes and protein-coding linc
standardized_tpm_expression_matrix_filtered_genes_file = pseudotissue_expression_dir + pseudotissue + '_normalized_expression_autosomes_protein_coding_and_linc.txt'
filter_expression_matrix(standardized_tpm_expression_matrix_file, standardized_tpm_expression_matrix_filtered_genes_file, gene_name_to_gene_info)



# Extract covariates (expression pcs)
num_samples = len(pseudotissue_ordered_sample_names)
if num_samples < 150:
	num_expression_pcs = 15
if num_samples >= 150 and num_samples < 250:
	num_expression_pcs = 30
elif num_samples >= 250 and num_samples < 350:
	num_expression_pcs = 45
elif num_samples >= 350 and num_samples < 1000:
	num_expression_pcs = 60
elif num_samples >= 1000:
	num_expression_pcs = 80
else:
	print('asssumption erororor')
	pdb.set_trace()
print(str(num_expression_pcs) + ' expression pcs selected because there are ' + str(num_samples) + ' samples')


# Generate expression pcs
standardized_tpm_expression_pcs_file = pseudotissue_expression_dir + pseudotissue + '_normalized_expression_pcs.txt'
generate_expression_pcs(standardized_tpm_expression_matrix_file, standardized_tpm_expression_pcs_file, num_expression_pcs)


