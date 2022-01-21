import numpy as np 
import os
import sys
import pdb





def extract_expression_pcs(pseudotissue_expr_pc_file, pseudotissue_sample_names):
	raw_expr_pc = np.loadtxt(pseudotissue_expr_pc_file,dtype=str, delimiter='\t')
	expr_pc = np.transpose(raw_expr_pc[1:,1:])
	expr_pc_sample_names = raw_expr_pc[1:,0]
	expr_pc_identifiers = raw_expr_pc[0,1:]

	if np.array_equal(np.asarray(expr_pc_sample_names), pseudotissue_sample_names) == False:
		print('assumption eororor')
		pdb.set_trace()
	new_expr_pc_identifiers = []
	for expr_pc_identifier in expr_pc_identifiers:
		pc_num = int(expr_pc_identifier.split('C')[1]) + 1
		new_identifier = 'expr_pc' + str(pc_num)
		new_expr_pc_identifiers.append(new_identifier)
	return expr_pc, np.asarray(new_expr_pc_identifiers)


def extract_list_of_covariates(gtex_covariate_dir, composit_tissues):
	num_tissues = len(composit_tissues)
	covariate_dicti = {}
	for composit_tissue in composit_tissues:
		covariate_file = gtex_covariate_dir + composit_tissue + '.v8.covariates.txt'
		covariate_data = np.loadtxt(covariate_file,dtype=str,delimiter='\t')
		covariate_names = covariate_data[1:,0]
		for covariate_name in covariate_names:
			if covariate_name.startswith('InferredCov'):
				continue
			if covariate_name not in covariate_dicti:
				covariate_dicti[covariate_name] = 0
			covariate_dicti[covariate_name] = covariate_dicti[covariate_name] + 1
	keys = np.sort([*covariate_dicti])
	valid_covariates = []
	for key in keys:
		if covariate_dicti[key] == num_tissues:
			valid_covariates.append(key)
	return np.sort(valid_covariates)

def extract_covariates_from_gtex_covariate_files(composit_tissues, gtex_covariate_dir, pseudotissue_sample_names):
	num_samples = len(pseudotissue_sample_names)
	# Create mapping from sample_name to matrix column index
	sample_mapping = {}
	for i, pseudotissue_sample_name in enumerate(pseudotissue_sample_names):
		sample_mapping[pseudotissue_sample_name] = i

	# Now extract list of all covariates (other than peer factors) found in all composit tissues
	covariates = extract_list_of_covariates(gtex_covariate_dir, composit_tissues)
	num_covariates = len(covariates)
	# Create mapping from covariate_name to matrix row index
	covariate_mapping = {}
	for i, covariate_name in enumerate(covariates):
		covariate_mapping[covariate_name] = i

	# Initialize covariate matrix
	cov_mat = np.zeros((num_covariates,num_samples)) - 1000

	for composit_tissue in composit_tissues:
		covariate_file = gtex_covariate_dir + composit_tissue + '.v8.covariates.txt'
		covariate_data = np.loadtxt(covariate_file,dtype=str,delimiter='\t')
		covariate_names = covariate_data[1:,0]
		indi_names = covariate_data[0,1:]
		cov = covariate_data[1:,1:]
		composit_tissue_sample_names = []
		for indi_name in indi_names:
			composit_tissue_sample_names.append(composit_tissue + ':' + indi_name)
		composit_tissue_sample_names = np.asarray(composit_tissue_sample_names)
		for ii, covariate_name in enumerate(covariate_names):
			for jj, composit_tissue_sample_name in enumerate(composit_tissue_sample_names):
				if covariate_name not in covariate_mapping:
					continue
				if composit_tissue_sample_name not in sample_mapping:
					continue

				row_index = covariate_mapping[covariate_name]
				col_index = sample_mapping[composit_tissue_sample_name]
				if cov_mat[row_index, col_index] != -1000.0:
					print('assumption eororor')
					pdb.set_trace()
				cov_mat[row_index, col_index] = float(cov[ii,jj])
	if np.sum(cov_mat == -1000.0) != 0:
		print('assumption erroro')
		pdb.set_trace()
	return cov_mat.astype(str), covariates	

### Command line args
pseudotissue = sys.argv[1]
composit_tissue_string = sys.argv[2]
pseudotissue_sample_names_file = sys.argv[3]
pseudotissue_expr_pc_file = sys.argv[4]
gtex_covariate_dir = sys.argv[5]
covariate_output_file = sys.argv[6]

# Extract list of composit tissues
composit_tissues = np.asarray(composit_tissue_string.split(','))

# Load in sample names
pseudotissue_sample_names = np.loadtxt(pseudotissue_sample_names_file,dtype=str, delimiter='\t')

# Extract expression pcs
expr_pc_mat, expr_pc_identifiers = extract_expression_pcs(pseudotissue_expr_pc_file, pseudotissue_sample_names)


# Extract other covariates
cov_mat, cov_identifiers = extract_covariates_from_gtex_covariate_files(composit_tissues, gtex_covariate_dir, pseudotissue_sample_names)


if cov_mat.shape[1] != expr_pc_mat.shape[1]:
	print('assumption eroror')
	pdb.set_trace()



# Print to output file
t = open(covariate_output_file,'w')
t.write('id\t' + '\t'.join(pseudotissue_sample_names) + '\n')

for index, expr_pc_identifier in enumerate(expr_pc_identifiers):
	arr = expr_pc_mat[index,:]
	if np.var(arr.astype(float)) == 0:
		print('assumption eroror')
		pdb.set_trace()
	t.write(expr_pc_identifier + '\t' + '\t'.join(expr_pc_mat[index,:]) + '\n')

for index, cov_identifier in enumerate(cov_identifiers):
	arr = cov_mat[index,:]
	if np.var(arr.astype(float)) == 0:
		print('assumption eroror')
		pdb.set_trace()
	t.write(cov_identifier + '\t' + '\t'.join(cov_mat[index,:]) + '\n')

t.close()



