import numpy as np 
import os
import sys
import pdb
import math
import scipy.stats
import coloc



def get_coloc_object(study_name, snp_names, sample_size, beta_vec, var_beta_vec):
	dicti = {}
	dicti['study_name'] = study_name
	dicti['snp_names'] = snp_names
	dicti['sample_size'] = sample_size
	dicti['beta'] = beta_vec
	dicti['var_beta'] = var_beta_vec
	return dicti


def write_mat_to_output_file(mat, row_names, col_names, output_file):
	t = open(output_file,'w')
	t.write('\t'.join(col_names) + '\n')
	for row_index in range(len(row_names)):
		t.write(row_names[row_index] + '\t' + '\t'.join(mat[row_index,:]) + '\n')
	t.close()

def run_coloc_for_single_gene_by_meta_analyzing_betas(gene_name, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name):
	# Coloc thresholds
	coloc_thresholds = np.asarray([.5, .7, .9, 95, .99])

	# Load in coloc data
	beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
	std_err_df = np.loadtxt(std_err_file, dtype=str,delimiter='\t')
	n_df = np.loadtxt(n_file, dtype=str,delimiter='\t')

	# Get studies observed for this gene
	observed_studies = beta_df[1:, 0]
	num_studies = len(observed_studies)
	eqtl_studies = observed_studies[:-1]

	# Get snp names
	snp_names = beta_df[0,1:]

	# Get object for trait for coloc
	trait_coloc_object = get_coloc_object(observed_studies[-1], snp_names, float(n_df[-1,1]), beta_df[-1,1:].astype(float), np.square(std_err_df[-1, 1:].astype(float)))


	# Make output matrix to pph values across tissues
	pph_mat = np.zeros((len(tissue_names), 5))
	pph_mat[:,0] = 1.0

	# Keep track of log bayes output matrix
	log_bayes_mat = np.reshape(['NULL']*(len(tissue_names)*5), (len(tissue_names), 5)).astype('U16')
	log_bayes_mat[:,-1] = str(len(snp_names))

	# Keep track of snp pph4s
	snp_pph4_mat = np.zeros((len(tissue_names), len(snp_names)))


	# Initilize output
	predicted_effects_list = []
	coloc_at_threshold_arr = []
	for coloc_threshold in coloc_thresholds:
		coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
		predicted_effects_list.append(coloc_thresh_mat)
		coloc_at_threshold_arr.append(False)

	trait_approx_lbf = coloc.get_approx_log_bf_estimates(trait_coloc_object)
	trait_beta = trait_coloc_object['beta']
	trait_var_beta = trait_coloc_object['var_beta']
	eqtl_approx_lbfs_arr = []
	observed_pph4s = []
	eqtl_beta_arr = []
	eqtl_var_beta_arr = []

	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(eqtl_studies):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]

		# Get object for eqtl study for coloc
		eqtl_coloc_object = get_coloc_object(observed_studies[eqtl_study_num], snp_names, float(n_df[(eqtl_study_num+1),1]), beta_df[(eqtl_study_num+1),1:].astype(float), np.square(std_err_df[(eqtl_study_num+1), 1:].astype(float)))

		# Estimate log bayes factors for each study
		eqtl_approx_lbf = coloc.get_approx_log_bf_estimates(eqtl_coloc_object)
		eqtl_approx_lbfs_arr.append(eqtl_approx_lbf)

		# Get beta and var beta
		eqtl_beta_arr.append(eqtl_coloc_object['beta'])
		eqtl_var_beta_arr.append(eqtl_coloc_object['var_beta'])
		
		# Esimate SNP PPH4 from log bayes factors
		snp_pph4 = coloc.get_snp_pph4_from_log_bf(eqtl_approx_lbf, trait_approx_lbf)
		snp_pph4[snp_pph4 < .01] = 0.0
		snp_pph4_mat[global_tissue_position,:] = snp_pph4

		# Get log bayes sums (basically summary stats for the gene
		lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = coloc.get_log_bayes_sums(eqtl_approx_lbf, trait_approx_lbf)

		# Keep track log sum bayes for this gene in each tissue
		log_bayes_mat[global_tissue_position,:-1] = np.asarray([lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4]).astype(str)

		# Keep track of Coloc probabilities
		pph_vec = coloc.run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)
		pph_mat[global_tissue_position,:] = pph_vec
		observed_pph4s.append(pph_vec[-1])

		# Compute predicted effects at each threshold
		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if pph_vec[4] > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4*pph_vec[4]*trait_coloc_object['beta']

	eqtl_approx_lbfs_arr = np.asarray(eqtl_approx_lbfs_arr)
	observed_pph4s = np.asarray(observed_pph4s)
	eqtl_beta_arr = np.asarray(eqtl_beta_arr)
	eqtl_var_beta_arr = np.asarray(eqtl_var_beta_arr)

	# Run competitive coloc for cases where gene has multiple tissues measured and at least one tissue has modest evidence for colocalization
	if len(eqtl_studies) > 1 and sum(pph_mat[:,-1] > .5) > 0:
		competitive_pph, study_posterior_probabilities, optimization_warning = coloc.competitive_coloc_study_posterior_probability_estimation_via_meta_analysis(trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr, observed_pph4s)
		if optimization_warning != 0:
			print('Scipy optimization did not converge successfully')
			print(gene_name)
		print('#######################')
		print(gene_name)
		print('max individual PPH4: ' + str(max(observed_pph4s)) + ' / competitive PPH4: ' + str(competitive_pph[-1]))
		for tissue_num, tissue_name in enumerate(eqtl_studies):
			print(str(observed_pph4s[tissue_num]) + '\t' + str(study_posterior_probabilities[tissue_num]))
		print(np.corrcoef(observed_pph4s, study_posterior_probabilities))


def run_coloc_for_single_gene_by_averaging_abfs(gene_name, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name):
	# Coloc thresholds
	coloc_thresholds = np.asarray([.1, .3, .5, .7, .9])

	# Load in coloc data
	beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
	std_err_df = np.loadtxt(std_err_file, dtype=str,delimiter='\t')
	n_df = np.loadtxt(n_file, dtype=str,delimiter='\t')

	# Get studies observed for this gene
	observed_studies = beta_df[1:, 0]
	num_studies = len(observed_studies)
	eqtl_studies = observed_studies[:-1]

	# Get snp names
	snp_names = beta_df[0,1:]

	# Get object for trait for coloc
	trait_coloc_object = get_coloc_object(observed_studies[-1], snp_names, float(n_df[-1,1]), beta_df[-1,1:].astype(float), np.square(std_err_df[-1, 1:].astype(float)))


	# Make output matrix to pph values across tissues
	pph_mat = np.zeros((len(tissue_names), 6))
	pph_mat[:,0] = 1.0

	pph_mat_v2 = np.zeros((len(tissue_names), 6))
	pph_mat_v2[:,0] = 1.0

	# Keep track of log bayes output matrix
	log_bayes_mat = np.reshape(['NULL']*(len(tissue_names)*5), (len(tissue_names), 5)).astype('U16')
	log_bayes_mat[:,-1] = str(len(snp_names))

	# Keep track of snp pph4s
	snp_pph4_mat = np.zeros((len(tissue_names), len(snp_names)))


	trait_approx_lbf = coloc.get_approx_log_bf_estimates(trait_coloc_object)
	eqtl_approx_lbfs_arr = []
	observed_pph4s = []

	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(eqtl_studies):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]

		# Get object for eqtl study for coloc
		eqtl_coloc_object = get_coloc_object(observed_studies[eqtl_study_num], snp_names, float(n_df[(eqtl_study_num+1),1]), beta_df[(eqtl_study_num+1),1:].astype(float), np.square(std_err_df[(eqtl_study_num+1), 1:].astype(float)))

		# Estimate log bayes factors for each study
		eqtl_approx_lbf = coloc.get_approx_log_bf_estimates(eqtl_coloc_object)
		eqtl_approx_lbfs_arr.append(eqtl_approx_lbf)
		
		# Esimate SNP PPH4 from log bayes factors
		snp_pph4 = coloc.get_snp_pph4_from_log_bf(eqtl_approx_lbf, trait_approx_lbf)
		snp_pph4[snp_pph4 < .001] = 0.0
		snp_pph4_mat[global_tissue_position,:] = snp_pph4

		# Get log bayes sums (basically summary stats for the gene
		lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = coloc.get_log_bayes_sums(eqtl_approx_lbf, trait_approx_lbf)

		# Keep track of Coloc probabilities
		pph_vec = coloc.run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)
		pph_mat[global_tissue_position,:-1] = pph_vec
		pph_mat_v2[global_tissue_position,:-1] = pph_vec
		observed_pph4s.append(pph_vec[-1])


	eqtl_approx_lbfs_arr = np.asarray(eqtl_approx_lbfs_arr)
	observed_pph4s = np.asarray(observed_pph4s)
	
	# Causal coloc probs
	causal_coloc_prob = coloc.get_causal_coloc_prob_from_pph4_vec(observed_pph4s)
	causal_coloc_prob_v2 = coloc.get_causal_coloc_prob_from_pph4_vec_v2(observed_pph4s)

	## GET PREDICTED EFFECTs FOR CAUSAL COLOC
	# Initilize output
	predicted_effects_list = []
	coloc_at_threshold_arr = []
	for coloc_threshold in coloc_thresholds:
		coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
		predicted_effects_list.append(coloc_thresh_mat)
		coloc_at_threshold_arr.append(False)
	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(eqtl_studies):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]
		pph_mat[global_tissue_position, -1] = causal_coloc_prob[eqtl_study_num]
		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if (observed_pph4s[eqtl_study_num]*causal_coloc_prob[eqtl_study_num]) > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4_mat[global_tissue_position,:]*trait_coloc_object['beta']
	# Save predicted effects mat
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		if coloc_at_threshold_arr[threshold_iter]:
			predicted_effect_size_file = coloc_output_dir + trait_name + '_' + gene_name + '_causal_coloc_' + str(coloc_threshold) + '_predicted_effect_sizes.txt'
			write_mat_to_output_file(predicted_effects_list[threshold_iter].astype(str), snp_names, tissue_names, predicted_effect_size_file)

	## GET PREDICTED EFFECTs FOR CAUSAL COLOC v2
	# Initilize output
	predicted_effects_list = []
	coloc_at_threshold_arr = []
	for coloc_threshold in coloc_thresholds:
		coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
		predicted_effects_list.append(coloc_thresh_mat)
		coloc_at_threshold_arr.append(False)
	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(eqtl_studies):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]
		pph_mat_v2[global_tissue_position, -1] = causal_coloc_prob_v2[eqtl_study_num]
		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if (observed_pph4s[eqtl_study_num]*causal_coloc_prob_v2[eqtl_study_num]) > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4_mat[global_tissue_position,:]*trait_coloc_object['beta']
	# Save predicted effects mat
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		if coloc_at_threshold_arr[threshold_iter]:
			predicted_effect_size_file = coloc_output_dir + trait_name + '_' + gene_name + '_causal_v2_coloc_' + str(coloc_threshold) + '_predicted_effect_sizes.txt'
			write_mat_to_output_file(predicted_effects_list[threshold_iter].astype(str), snp_names, tissue_names, predicted_effect_size_file)


	## GET PREDICTED EFFECTs FOR STANDARD COLOC
	# Initilize output
	predicted_effects_list = []
	coloc_at_threshold_arr = []
	for coloc_threshold in coloc_thresholds:
		coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
		predicted_effects_list.append(coloc_thresh_mat)
		coloc_at_threshold_arr.append(False)
	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(eqtl_studies):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]
		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if (observed_pph4s[eqtl_study_num]) > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4_mat[global_tissue_position,:]*trait_coloc_object['beta']
	# Save predicted effects mat
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		if coloc_at_threshold_arr[threshold_iter]:
			predicted_effect_size_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_' + str(coloc_threshold) + '_predicted_effect_sizes.txt'
			write_mat_to_output_file(predicted_effects_list[threshold_iter].astype(str), snp_names, tissue_names, predicted_effect_size_file)


	# Save snp pph4 mat to output file
	snp_pph4_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_snp_pph4.txt'
	write_mat_to_output_file(snp_pph4_mat.astype(str), tissue_names, snp_names, snp_pph4_output_file)

	# Save pph mat
	pph_mat_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
	write_mat_to_output_file(pph_mat.astype(str), tissue_names, np.asarray(["PPH0", "PPH1", "PPH2", "PPH3", "PPH4", "PPHC"]), pph_mat_output_file)

	# Save pph mat
	pph_mat_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_v2_coloc_posterior_probabilities.txt'
	write_mat_to_output_file(pph_mat_v2.astype(str), tissue_names, np.asarray(["PPH0", "PPH1", "PPH2", "PPH3", "PPH4", "PPHC"]), pph_mat_output_file)


gene_file = sys.argv[1]
trait_name = sys.argv[2]
gtex_tissue_file = sys.argv[3]
coloc_output_dir = sys.argv[4]
job_number = float(sys.argv[5])
total_jobs = float(sys.argv[6])



# Load in tissue names
tissue_df = np.loadtxt(gtex_tissue_file,dtype=str,delimiter='\t')
tissue_names = tissue_df[1:,0]
# Create mapping from tissue name to position
tissue_name_to_position = {}
for index, tissue_name in enumerate(tissue_names):
	tissue_name_to_position[tissue_name] = index


# Load in gene data frame
gene_df = np.loadtxt(gene_file, dtype=str,delimiter='\t')
gene_df_header = gene_df[0,:]
gene_df = gene_df[1:,:]
# Get total number of genes
num_genes = gene_df.shape[0]


# For parallelization purposes, determine which genes to test in this thread
tasks_per_job = np.floor(num_genes/total_jobs) + 1
start_task = int(np.floor(job_number*tasks_per_job + 1) - 1)
end_task = int(np.floor((job_number + 1)*tasks_per_job) -1)
if end_task > (num_genes - 1):
	end_task = num_genes -1




# Loop through all genes and run colocalization for that gene
for gene_num in range(start_task, (end_task+1)):
	#print(gene_num)
	gene_name = gene_df[gene_num, 0]
	n_file = gene_df[gene_num, 4]
	beta_file = gene_df[gene_num, 5]
	std_err_file = gene_df[gene_num, 6]

	# Run coloc for a single gene
	run_coloc_for_single_gene_by_averaging_abfs(gene_name, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name)
	#run_coloc_for_single_gene_by_meta_analyzing_betas(gene_name, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name)

