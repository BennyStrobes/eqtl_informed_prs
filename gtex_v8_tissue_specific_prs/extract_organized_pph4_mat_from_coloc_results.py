import numpy as np 
import os
import sys
import pdb



def load_in_pph4s_across_genes(gene_names, trait_name, coloc_output_dir):
	pph4_mat = []
	counter = 0
	for gene_name in gene_names:
		counter = counter + 1
		gene_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		data = np.loadtxt(gene_pph4_file,dtype=str,delimiter='\t', skiprows=1)
		pph4 = data[:,5].astype(float)
		pph4_mat.append(pph4)
	return np.asarray(pph4_mat)

def load_in_pphcs_across_genes_v1(gene_names, trait_name, coloc_output_dir):
	pph4_mat = []
	counter = 0
	for gene_name in gene_names:
		counter = counter + 1
		gene_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		data = np.loadtxt(gene_pph4_file,dtype=str,delimiter='\t', skiprows=1)
		pph4 = data[:,6].astype(float)
		pph4_mat.append(pph4)
	return np.asarray(pph4_mat)

def load_in_pphcs_across_genes_v2(gene_names, trait_name, coloc_output_dir):
	pph4_mat = []
	counter = 0
	for gene_name in gene_names:
		counter = counter + 1
		gene_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_v2_coloc_posterior_probabilities.txt'
		data = np.loadtxt(gene_pph4_file,dtype=str,delimiter='\t', skiprows=1)
		pph4 = data[:,6].astype(float)
		pph4_mat.append(pph4)
	return np.asarray(pph4_mat)


def load_in_betas_across_genes(gene_names, trait_name, coloc_output_dir, gene_df, col_index):
	betas = []
	for gene_index, gene_name in enumerate(gene_names):
		print(gene_index)
		beta_file = gene_df[gene_index,col_index]
		beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
		snp_names = beta_df[0,1:]
		beta_vec = beta_df[-1,1:].astype(float)
		snp_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_snp_pph4.txt'
		snp_pph4 = np.loadtxt(snp_pph4_file,dtype=str,delimiter='\t',skiprows=1)
		snp_pph4 = snp_pph4[:,1:].astype(float)
		num_tiss = snp_pph4.shape[0]
		if num_tiss != 23:
			print('assumption erororo')
			pdb.set_trace()
		top_beta = np.zeros(23)
		for tissue_num in range(num_tiss):
			top_beta[tissue_num] = beta_vec[np.argmax(snp_pph4[tissue_num,:])]
		betas.append(top_beta)
	return np.asarray(betas)

def get_number_of_times_each_variant_is_loaded(gene_names, trait_name, coloc_output_dir, gene_df, tissue_names, tissue_name, pph4_threshold):
	tissue_index = np.where(tissue_names==tissue_name)[0][0]
	dicti = {}
	for gene_index, gene_name in enumerate(gene_names):
		print(gene_index)
		# gene pph4
		gene_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		data = np.loadtxt(gene_pph4_file,dtype=str,delimiter='\t', skiprows=1)
		pph4 = data[:,6].astype(float)
		tissue_pph4 = pph4[tissue_index]
		if tissue_pph4 < pph4_threshold:
			continue
		# Betas
		beta_file = gene_df[gene_index,5]
		beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
		snp_names = beta_df[0,1:]
		beta_vec = beta_df[-1,1:].astype(float)
		# Snp pph4
		snp_pph4_file = coloc_output_dir + trait_name + '_' + gene_name + '_snp_pph4.txt'
		snp_pph4 = np.loadtxt(snp_pph4_file,dtype=str,delimiter='\t',skiprows=1)
		tissue_snp_pph4 = snp_pph4[tissue_index,1:].astype(float)
		top_snp = np.argmax(tissue_snp_pph4)
		#snp_indices = tissue_snp_pph4>.01

		# Select snps loaded on this gene
		snp_name_on_gene = snp_names[top_snp]
		snp_beta_on_gene = beta_vec[top_snp]

		if snp_name_on_gene not in dicti:
			dicti[snp_name_on_gene] = (1, snp_beta_on_gene)
		else:
			old_tuple = dicti[snp_name_on_gene]
			dicti[snp_name_on_gene] = (old_tuple[0] + 1, snp_beta_on_gene)
	
	keys = [*dicti]
	t = open(coloc_output_dir + trait_name + '_snp_repeats_in_' + tissue_name + '.txt','w')
	t.write('snp_name\tN\tsnp_beta\n')
	for snp_name in keys:
		t.write(snp_name + '\t' + str(dicti[snp_name][0]) + '\t' +str(dicti[snp_name][1]) + '\n')
	t.close()

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

# For a single tissue, export number of times each variant is loaded on a gene with pph4 > threshold
tissue_name = 'Whole_Blood'
pph4_threshold = .5
get_number_of_times_each_variant_is_loaded(gene_names, trait_name, coloc_output_dir, gene_df, tissue_names, tissue_name, pph4_threshold)



'''
# PPCs
pph4_mat = load_in_pphcs_across_genes_v1(gene_names, trait_name, coloc_output_dir)
t = open(coloc_output_dir + trait_name + '_pphc_v1_matrix.txt','w')
t.write('gene_name\t' + '\t'.join(tissue_names) + '\n')
for gene_index, gene_name in enumerate(gene_names):
	t.write(gene_name + '\t' + '\t'.join(pph4_mat[gene_index,:].astype(str)) + '\n')
t.close()

# PPCs
pph4_mat = load_in_pphcs_across_genes_v2(gene_names, trait_name, coloc_output_dir)
t = open(coloc_output_dir + trait_name + '_pphc_v2_matrix.txt','w')
t.write('gene_name\t' + '\t'.join(tissue_names) + '\n')
for gene_index, gene_name in enumerate(gene_names):
	t.write(gene_name + '\t' + '\t'.join(pph4_mat[gene_index,:].astype(str)) + '\n')
t.close()



# PPH4s
pph4_mat = load_in_pph4s_across_genes(gene_names, trait_name, coloc_output_dir)
t = open(coloc_output_dir + trait_name + '_pph4_matrix.txt','w')
t.write('gene_name\t' + '\t'.join(tissue_names) + '\n')
for gene_index, gene_name in enumerate(gene_names):
	t.write(gene_name + '\t' + '\t'.join(pph4_mat[gene_index,:].astype(str)) + '\n')
t.close()



# Z scores
z_mat = load_in_betas_across_genes(gene_names, trait_name, coloc_output_dir, gene_df, 3)
t = open(coloc_output_dir + trait_name + '_top_z_matrix.txt','w')
t.write('gene_name\t' + '\t'.join(tissue_names) + '\n')
for gene_index, gene_name in enumerate(gene_names):
	t.write(gene_name + '\t' + '\t'.join(z_mat[gene_index,:].astype(str)) + '\n')
t.close()



# BETAS
beta_mat = load_in_betas_across_genes(gene_names, trait_name, coloc_output_dir, gene_df, 5)
t = open(coloc_output_dir + trait_name + '_top_beta_matrix.txt','w')
t.write('gene_name\t' + '\t'.join(tissue_names) + '\n')
for gene_index, gene_name in enumerate(gene_names):
	t.write(gene_name + '\t' + '\t'.join(beta_mat[gene_index,:].astype(str)) + '\n')
t.close()
'''

