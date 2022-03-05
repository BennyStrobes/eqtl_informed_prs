import numpy as np 
import os
import sys
import pdb
import coloc



def extract_lb_mat_data_for_single_tissue(tissue_name, tissue_position, gene_names, trait_name, coloc_output_dir):
	lb_mat = []
	num_snps_arr = []
	counter = 0
	for gene_name in gene_names:
		counter = counter + 1
		gene_lb_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_log_bayes_sums.txt'
		data = np.loadtxt(gene_lb_file,dtype=str,delimiter='\t', skiprows=1)
		if data[tissue_position,0] != tissue_name:
			print('assumption erooror')
			pdb.set_trace()
		if data[tissue_position,1] == 'NULL':
			continue
		line_lb = data[tissue_position,1:5].astype(float)
		num_snps = data[tissue_position,5].astype(float)
		lb_mat.append(line_lb)
		num_snps_arr.append(num_snps)
	return np.asarray(lb_mat), np.asarray(num_snps_arr)

def update_gene_pph_mat_given_fixed_coloc_prior(lb_mat, coloc_prior):
	num_genes = lb_mat.shape[0]
	pph_mat = []
	for gene_num in range(num_genes):
		pph_vec = coloc.run_coloc_with_precomputed_log_bayes_sums(lb_mat[gene_num,0], lb_mat[gene_num,1], lb_mat[gene_num,2], lb_mat[gene_num,3], p0=coloc_prior[0], p1=coloc_prior[1], p2=coloc_prior[2], p12=coloc_prior[3])
		pph_mat.append(pph_vec)
	return np.asarray(pph_mat)

def update_coloc_prior_given_fixed_pph_mat(pph_mat, num_snps, pseudocounts):
	un_normalized_coloc_prior = np.copy(pseudocounts)

	# Get number of genes for this tissue
	num_genes = pph_mat.shape[0]
	if num_genes != len(num_snps):
		print('assumption eroror')
		pdb.set_trace()

	# Loop through genes
	for gene_num in range(num_genes):
		gene_pph_vec = pph_mat[gene_num,:]
		gene_num_snps = num_snps[gene_num]

		# Count up expected number of snps associated with neither trait 1 or trait 2
		un_normalized_coloc_prior[0] = un_normalized_coloc_prior[0] + (gene_num_snps*gene_pph_vec[0]) + ((gene_num_snps-1.0)*gene_pph_vec[1]) + ((gene_num_snps-1.0)*gene_pph_vec[2]) + ((gene_num_snps-2.0)*gene_pph_vec[3]) + ((gene_num_snps-1.0)*gene_pph_vec[4])
		# Count up expected number of snps ONLY associated with trait 1
		un_normalized_coloc_prior[1] = un_normalized_coloc_prior[1] + gene_pph_vec[1] + gene_pph_vec[3]
		# Count up expected number of snps ONLY associated with trait 2
		un_normalized_coloc_prior[2] = un_normalized_coloc_prior[2] + gene_pph_vec[2] + gene_pph_vec[3]
		# Count up expected number of snps associated with BOTH trait 1 and trait 2
		un_normalized_coloc_prior[3] = un_normalized_coloc_prior[3] + gene_pph_vec[4]
	
	return un_normalized_coloc_prior/np.sum(un_normalized_coloc_prior)


def estimate_coloc_priors_via_iterative_algorithm(lb_mat, num_snps, p1_alpha=1e-4, p2_alpha=1e-4, p12_alpha=1e-5, scale=100.0, max_iterations=50):
	p0_alpha = 1.0 - p1_alpha - p2_alpha - p12_alpha
	coloc_prior_init = np.asarray([p0_alpha, p1_alpha, p2_alpha, p12_alpha])
	pseudocounts = coloc_prior_init*scale
	coloc_prior = np.copy(coloc_prior_init)

	# Begin iterative algorithm
	for iterative_iteration in range(max_iterations):
		# First update gene PPHs given fixed coloc_prior
		pph_mat = update_gene_pph_mat_given_fixed_coloc_prior(lb_mat, coloc_prior)
		# Then update coloc prior given fixed gene PPHs
		coloc_prior = update_coloc_prior_given_fixed_pph_mat(pph_mat, num_snps, pseudocounts)

	return coloc_prior


def write_mat_to_output_file(mat, row_names, col_names, output_file):
	t = open(output_file,'w')
	t.write('\t'.join(col_names) + '\n')
	for row_index in range(len(row_names)):
		t.write(row_names[row_index] + '\t' + '\t'.join(mat[row_index,:]) + '\n')
	t.close()


def run_causal_coloc_for_single_gene(tissue_names, tissue_name_to_position, gene_name, causal_prior, beta_file, snp_pph4_file, gene_pph4_file, coloc_output_dir, trait_name):
	coloc_thresholds = np.asarray([.1, .3, .5, .7, .9])
	# Make output matrix to pph values across tissues
	pph_mat = np.zeros((len(tissue_names), 6))
	pph_mat[:,0] = 1.0

	pp_data = np.loadtxt(gene_pph4_file,dtype=str,delimiter='\t', skiprows=1)

	pph4_vec = pp_data[:,5].astype(float)

	gene_posteriors = coloc.get_causal_coloc_prob_from_pph4_vec_v2_with_prior(pph4_vec, causal_prior)

	pp_data[:,6] = gene_posteriors.astype(str)

	pph_mat = pp_data[:,1:].astype(float)


	# MM coloc
	boolean = False
	coloc_at_threshold_arr = []
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		coloc_at_threshold_arr.append(False)
	coloc_at_threshold_arr = np.asarray(coloc_at_threshold_arr)

	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(tissue_names):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]


		if pph4_vec[global_tissue_position]*gene_posteriors[global_tissue_position] > np.min(coloc_thresholds) and boolean == False:
			boolean = True
			beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
			snp_names = beta_df[0,1:]
			beta_vec = beta_df[-1,1:].astype(float)

			snp_pph4 = np.loadtxt(snp_pph4_file,dtype=str,delimiter='\t',skiprows=1)
			snp_pph4 = snp_pph4[:,1:].astype(float)
			# Initilize output
			predicted_effects_list = []
			for coloc_threshold in coloc_thresholds:
				coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
				predicted_effects_list.append(coloc_thresh_mat)



		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if pph4_vec[global_tissue_position]*gene_posteriors[global_tissue_position] > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4[global_tissue_position,:]*pph4_vec[global_tissue_position]*gene_posteriors[global_tissue_position]*beta_vec

	# Save predicted effects mat
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		if coloc_at_threshold_arr[threshold_iter]:
			predicted_effect_size_file = coloc_output_dir + trait_name + '_' + gene_name + '_mm_coloc_' + str(coloc_threshold) + '_predicted_effect_sizes.txt'
			write_mat_to_output_file(predicted_effects_list[threshold_iter].astype(str), snp_names, tissue_names, predicted_effect_size_file)


	# MM coloc non causal thresh
	boolean = False
	coloc_at_threshold_arr = []
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		coloc_at_threshold_arr.append(False)
	coloc_at_threshold_arr = np.asarray(coloc_at_threshold_arr)

	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(tissue_names):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]


		if pph4_vec[global_tissue_position] > np.min(coloc_thresholds) and boolean == False:
			boolean = True
			beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
			snp_names = beta_df[0,1:]
			beta_vec = beta_df[-1,1:].astype(float)

			snp_pph4 = np.loadtxt(snp_pph4_file,dtype=str,delimiter='\t',skiprows=1)
			snp_pph4 = snp_pph4[:,1:].astype(float)
			# Initilize output
			predicted_effects_list = []
			for coloc_threshold in coloc_thresholds:
				coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
				predicted_effects_list.append(coloc_thresh_mat)

		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if pph4_vec[global_tissue_position] > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4[global_tissue_position,:]*pph4_vec[global_tissue_position]*gene_posteriors[global_tissue_position]*beta_vec

	# Save predicted effects mat
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		if coloc_at_threshold_arr[threshold_iter]:
			predicted_effect_size_file = coloc_output_dir + trait_name + '_' + gene_name + '_mm_coloc_nc_thresh_' + str(coloc_threshold) + '_predicted_effect_sizes.txt'
			write_mat_to_output_file(predicted_effects_list[threshold_iter].astype(str), snp_names, tissue_names, predicted_effect_size_file)


	# Save pph mat
	pph_mat_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_mm_coloc_posterior_probabilities.txt'
	write_mat_to_output_file(pph_mat.astype(str), tissue_names, np.asarray(["PPH0", "PPH1", "PPH2", "PPH3", "PPH4", "PPHC"]), pph_mat_output_file)




def load_in_pph4s_across_genes(gene_names, trait_name, coloc_output_dir):
	pph4_mat = []
	counter = 0
	for gene_name in gene_names:
		counter = counter + 1
		gene_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		data = np.loadtxt(gene_pph4_file,dtype=str,delimiter='\t', skiprows=1)
		pph4 = data[:,5].astype(float)
		if np.max(pph4) < .5:
			continue
		if np.sum(pph4==0) > 0:
			continue
		pph4_mat.append(pph4)
	return np.asarray(pph4_mat)

def get_causal_posteriors(pph4_mat, causal_prior):
	num_genes = pph4_mat.shape[0]
	causal_posteriors_arr = []
	for gene_num in range(num_genes):
		gene_posteriors = coloc.get_causal_coloc_prob_from_pph4_vec_v2_with_prior(pph4_mat[gene_num,:], causal_prior)
		causal_posteriors_arr.append(gene_posteriors)
	return np.asarray(causal_posteriors_arr)

def run_iterative_algorithm_to_get_causal_prior(pph4_mat, pseudocount=1e-9, max_iters=100):
	num_dim = pph4_mat.shape[1]
	# Initialize causal prior
	causal_prior = np.zeros(num_dim) + 1.0/num_dim

	for itera in range(max_iters):
		causal_posteriors = get_causal_posteriors(pph4_mat, causal_prior)
		un_normalized = np.sum(causal_posteriors,axis=0) + pseudocount
		causal_prior = un_normalized/np.sum(un_normalized)
		print(causal_prior)
	return causal_prior


# Command line args
gene_file = sys.argv[1]
trait_name = sys.argv[2]
gtex_tissue_file = sys.argv[3]
coloc_output_dir = sys.argv[4]



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
gene_names = gene_df[:,0]

pph4_mat = load_in_pph4s_across_genes(gene_names, trait_name, coloc_output_dir)

causal_prior = run_iterative_algorithm_to_get_causal_prior(pph4_mat)



t = open(coloc_output_dir + trait_name + '_coloc_mm_causal_prior_probability.txt','w')
t.write('tissue_name\tprior_probability\n')
for index, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(causal_prior[index]) + '\n')
t.close()

# Loop through all genes and run colocalization for that gene
for gene_num in range(num_genes):
	gene_name = gene_df[gene_num, 0]
	print(gene_name)
	beta_file = gene_df[gene_num, 5]
	snp_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_snp_pph4.txt'
	gene_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
	run_causal_coloc_for_single_gene(tissue_names, tissue_name_to_position, gene_name, causal_prior, beta_file, snp_pph4_file, gene_pph4_file, coloc_output_dir, trait_name)
