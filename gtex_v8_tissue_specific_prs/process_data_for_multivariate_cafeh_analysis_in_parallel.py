import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import gzip
import pandas as pd


def get_gtex_tissue_names(gtex_tissue_file):
	f = open(gtex_tissue_file)
	head_count = 0
	tissues = []
	sample_sizes = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissues.append(data[0])
		sample_sizes.append(data[1])
	f.close()
	return np.asarray(tissues), np.asarray(sample_sizes)

def create_gwas_snp_id_to_beta_std_err_mapping_on_single_chromosome(sumstat_file, chrom_string, trait_sample_size):
	f = gzip.open(sumstat_file)
	dicti = {}
	error_counter = 0
	alt_dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip().decode('utf-8')
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Error checking
		if len(data) != 16:
			print('assumption eororoor')
			pdb.set_trace()
		# Skip snps not on desired chromosome
		line_string = data[1]
		if line_string != chrom_string:
			continue
		# Extract relevent fields (Effect allele is always last)
		snp_id = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4]
		alt_snp_id = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5]

		# Get allele frequency
		af = float(data[6])
		info = float(data[7])
		beta_bolt_lmm_inf = float(data[10])
		se_bolt_lmm_inf = float(data[11])
		chi_sq_bolt_lmm = float(data[14])

		if chi_sq_bolt_lmm < 0:
			error_counter = error_counter +1
			continue

		# Get maf from af
		if af > .5:
			maf = 1.0 - af
		else:
			maf = af
		if snp_id in dicti or alt_snp_id in alt_dicti:
			print('assumption oororor')
			pdb.set_trace()
		dicti[snp_id] = (beta_bolt_lmm_inf, se_bolt_lmm_inf, maf, trait_sample_size, info)
		alt_dicti[alt_snp_id] = (beta_bolt_lmm_inf, se_bolt_lmm_inf, maf, trait_sample_size, info)
	f.close()
	return dicti, alt_dicti


def get_genes_on_this_chromosome(cafeh_gene_list_file, chrom_string):
	gene_dicti = {}
	f = open(cafeh_gene_list_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_chrom_num = data[1]
		# remove genes not on this chromosome
		if line_chrom_num != chrom_string:
			continue
		gene_id = data[0]
		tissues = data[2].split(';')
		if gene_id in gene_dicti:
			print('assumptionototoe')
			pdb.set_trace()
		gene_dicti[gene_id] = tissues
	f.close()
	return gene_dicti

def extract_gene_level_gtex_data(gene_level_gtex_associations_file):
	f = open(gene_level_gtex_associations_file)
	head_count = 0
	snp_arr = []
	tissue_arr = []
	zscore_arr = []
	beta_arr = []
	std_err_arr = []
	snp_dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		# Quick error checking
		if len(data) != 6:
			continue
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[0]
		snp_id = data[2]
		beta = float(data[4])
		std_err = float(data[5])
		zscore=beta/std_err
		snp_arr.append(snp_id)
		zscore_arr.append(zscore)
		tissue_arr.append(tissue_name)
		beta_arr.append(beta)
		std_err_arr.append(std_err)

		if snp_id not in snp_dicti:
			snp_dicti[snp_id] = 0
		snp_dicti[snp_id] = snp_dicti[snp_id] + 1
	f.close()

	# Put in np arrays
	snp_arr = np.asarray(snp_arr)
	zscore_arr = np.asarray(zscore_arr)
	tissue_arr = np.asarray(tissue_arr)
	beta_arr = np.asarray(beta_arr)
	std_err_arr = np.asarray(std_err_arr)

	# Get list of snps (valid_snps) found in all tissues
	unique_tissues = np.unique(tissue_arr)
	num_unique_tissues = len(unique_tissues)
	unique_snps = np.asarray([*snp_dicti])

	valid_snps = {}
	for snp_name in unique_snps:
		if snp_dicti[snp_name] == num_unique_tissues:
			valid_snps[snp_name] = 1
	
	# Filter data to valid snps
	valid_indices = []
	for index, snp_name in enumerate(snp_arr):
		if snp_name in valid_snps:
			valid_indices.append(index)
	valid_indices = np.asarray(valid_indices)

	return snp_arr[valid_indices], zscore_arr[valid_indices], tissue_arr[valid_indices], beta_arr[valid_indices], std_err_arr[valid_indices]

def get_bivariate_zscores_matrix_and_sample_size_matrix(tissue_snps, tissue_zscores, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt, gtex_tissue, trait_name, gtex_tissue_sample_size):
	bivariate_zscores = []
	bivariate_betas = []
	bivariate_std_errs = []
	snp_names = []
	bivariate_sample_sizes = []
	for snp_index, snp_name in enumerate(tissue_snps):
		tissue_zscore = tissue_zscores[snp_index]
		if snp_name in gwas_snp_id_to_z_score:
			if gwas_snp_id_to_z_score[snp_name][2] < .001 or gwas_snp_id_to_z_score[snp_name][4] < .6:
				continue
			gwas_beta = gwas_snp_id_to_z_score[snp_name][0]
			gwas_std_err = gwas_snp_id_to_z_score[snp_name][1]
			gwas_zscore = gwas_beta/gwas_std_err
			bivariate_zscores.append((tissue_zscore, gwas_zscore))
			bivariate_betas.append(gwas_beta)
			bivariate_std_errs.append(gwas_std_err)
			snp_names.append(snp_name)
			bivariate_sample_sizes.append((gtex_tissue_sample_size, gwas_snp_id_to_z_score[snp_name][3]))
		elif snp_name in gwas_snp_id_to_z_score_alt:
			# ignore snp if maf < .1% or info score < .6
			if gwas_snp_id_to_z_score_alt[snp_name][2] < .001 or gwas_snp_id_to_z_score_alt[snp_name][4] < .6:
				continue
			gwas_beta = -1.0*gwas_snp_id_to_z_score_alt[snp_name][0]
			gwas_std_err = gwas_snp_id_to_z_score_alt[snp_name][1]
			gwas_zscore = gwas_beta/gwas_std_err			
			bivariate_zscores.append((tissue_zscore, gwas_zscore))
			bivariate_betas.append(gwas_beta)
			bivariate_std_errs.append(gwas_std_err)
			snp_names.append(snp_name)
			bivariate_sample_sizes.append((gtex_tissue_sample_size, gwas_snp_id_to_z_score_alt[snp_name][3]))
		else:
			continue
	bivariate_zscore_mat = np.asarray(bivariate_zscores)
	gwas_beta_arr = np.asarray(bivariate_betas)
	gwas_std_err_arr = np.asarray(bivariate_std_errs)
	snp_names = np.asarray(snp_names)
	bivariate_sample_size_mat = np.asarray(bivariate_sample_sizes)
	if len(snp_names) > 9:
		study_ids = np.asarray([gtex_tissue, trait_name])
		z_df = pd.DataFrame(np.transpose(bivariate_zscore_mat), index=study_ids, columns=snp_names)
		beta_df = pd.DataFrame(np.asmatrix(gwas_beta_arr), index=study_ids[1:], columns=snp_names)
		std_err_df = pd.DataFrame(np.asmatrix(gwas_std_err_arr), index=study_ids[1:], columns=snp_names)
		n_df = pd.DataFrame(np.transpose(bivariate_sample_size_mat.astype(int)), index=study_ids, columns=snp_names)
		boolers = True
	else:
		print('assumption error: Gene thrown out')
		z_df = 0
		n_df = 0
		beta_df = 0
		std_err_df = 0
		snp_names = 0
		boolers = False
	return z_df, beta_df, std_err_df, n_df, snp_names, boolers

def create_mapping_from_hg38_snp_id_to_hg19_snp_id(hg38_to_hg19_mapping_file):
	mapping = {}
	r_mapping = {}
	f = open(hg38_to_hg19_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		hg38_snp_id = data[0]
		hg19_snp_id = data[1]
		if hg19_snp_id == 'NA':
			continue
		mapping[hg38_snp_id] = hg19_snp_id
		r_mapping[hg19_snp_id] = hg38_snp_id
	f.close()
	return mapping, r_mapping

def mean_impute_missing_snps(re_ordered_geno):
	num_snps = re_ordered_geno.shape[1]
	new_arr = []
	for snp_num in range(num_snps):
		NA_indices = re_ordered_geno[:,snp_num]== 'NA'
		num_missing_indices = sum(NA_indices)
		if num_missing_indices == 0:
			new_arr.append(re_ordered_geno[:,snp_num].astype(float))
		else:
			observed_indices = re_ordered_geno[:,snp_num] != 'NA'
			meany = np.mean(re_ordered_geno[observed_indices,snp_num].astype(float))
			re_ordered_geno[NA_indices, snp_num] = str(meany)
			new_arr.append(re_ordered_geno[:,snp_num].astype(float))
	return np.transpose(np.asarray(new_arr))

def convert_genotype_file_to_hg19_and_reorder(genotype_file, hg38_to_hg19_snp_id, hg19_to_hg38_snp_id, gene_gtex_unique_snps):
	geno = np.loadtxt(genotype_file, dtype=str)
	donor_id = geno[0,1:]
	hg38_snps = geno[1:,0]
	geno = np.transpose(geno[1:,1:])

	mapping = {}
	hg19_snps = []
	for i,hg38_snp in enumerate(hg38_snps):
		hg38_snp = '_'.join(hg38_snp.split('_')[:5])
		if hg38_snp not in hg38_to_hg19_snp_id:
			hg19_snps.append('NA')
			continue
		hg19_snp = hg38_to_hg19_snp_id[hg38_snp]
		mapping[hg19_snp] = i
		hg19_snps.append(hg19_snp)
	hg19_snps = np.asarray(hg19_snps)
	if len(hg19_snps) != len(hg38_snps):
		print('assumption erororo')
		pdb.set_trace()

	indices = []
	for gtex_snp in gene_gtex_unique_snps:
		indices.append(mapping[gtex_snp])
	indices = np.asarray(indices)

	if np.array_equal(hg19_snps[indices],gene_gtex_unique_snps) == False:
		print('assumption error')
		pdb.set_trace()

	re_ordered_geno = geno[:, indices]
	g_df = pd.DataFrame(re_ordered_geno.astype(float), index=donor_id, columns=gene_gtex_unique_snps)

	return g_df

def filter_gtex_snps_to_those_that_have_gwas_support(gene_gtex_unique_snps, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt):
	cafeh_snps = []
	for snp_name in gene_gtex_unique_snps:
		if snp_name in gwas_snp_id_to_z_score:
			if gwas_snp_id_to_z_score[snp_name][2] >= .001 and gwas_snp_id_to_z_score[snp_name][4] >= .6:
				cafeh_snps.append(snp_name)
		elif snp_name in gwas_snp_id_to_z_score_alt:
			# ignore snp if maf < .1% or info score < .6
			if gwas_snp_id_to_z_score_alt[snp_name][2] >= .001 and gwas_snp_id_to_z_score_alt[snp_name][4] >= .6:
				cafeh_snps.append(snp_name)

	cafeh_snps = np.sort(np.asarray(cafeh_snps))

	if len(cafeh_snps) != len(np.unique(cafeh_snps)):
		print('assumption eroror')

	return cafeh_snps

def get_summary_stat_matrices(cafeh_snps, cafeh_tissues, gene_gtex_snps, gene_gtex_tissues, gene_gtex_zscores, gene_gtex_betas, gene_gtex_std_errs, trait_name, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt):
	num_snps = len(cafeh_snps)
	num_tissues = len(cafeh_tissues)

	cafeh_studies = np.append(cafeh_tissues, trait_name)

	# Create mapping from snp to row num
	snp_to_row_num = {}
	for index, snp_id in enumerate(cafeh_snps):
		snp_to_row_num[snp_id] = index

	# Create mapping from tissue to col num
	tissue_to_col_num = {}
	for index, tissue_id in enumerate(cafeh_tissues):
		tissue_to_col_num[tissue_id] = index

	z_score_mat = -500.0*np.ones((num_snps, (num_tissues+1)))
	n_mat = 320.0*np.ones((num_snps, (num_tissues+1)))

	for gtex_snp_index, gtex_snp_id in enumerate(gene_gtex_snps):
		gtex_tissue = gene_gtex_tissues[gtex_snp_index]
		gtex_zscore = gene_gtex_zscores[gtex_snp_index]

		if gtex_snp_id not in snp_to_row_num:
			continue
		if gtex_tissue not in tissue_to_col_num:
			continue

		z_score_mat[snp_to_row_num[gtex_snp_id], tissue_to_col_num[gtex_tissue]] = gtex_zscore
	
	beta_arr = []
	std_err_arr = []
	for snp_index, snp_id in enumerate(cafeh_snps):
		if snp_id in gwas_snp_id_to_z_score:
			if gwas_snp_id_to_z_score[snp_id][2] < .001 or gwas_snp_id_to_z_score[snp_id][4] < .6:
				continue
			gwas_beta = gwas_snp_id_to_z_score[snp_id][0]
			gwas_std_err = gwas_snp_id_to_z_score[snp_id][1]
			gwas_zscore = gwas_beta/gwas_std_err
			gwas_sample_size = gwas_snp_id_to_z_score[snp_id][3]
		elif snp_id in gwas_snp_id_to_z_score_alt:
			if gwas_snp_id_to_z_score_alt[snp_id][2] < .001 or gwas_snp_id_to_z_score_alt[snp_id][4] < .6:
				continue
			gwas_beta = -1.0*gwas_snp_id_to_z_score_alt[snp_id][0]
			gwas_std_err = gwas_snp_id_to_z_score_alt[snp_id][1]
			gwas_zscore = gwas_beta/gwas_std_err
			gwas_sample_size = gwas_snp_id_to_z_score_alt[snp_id][3]
		else:
			print('asssumption eroror: already checked to make sure snp in gwas')
			pdb.set_trace()

		z_score_mat[snp_index, -1] = gwas_zscore
		n_mat[snp_index, -1] = gwas_sample_size
		beta_arr.append(gwas_beta)
		std_err_arr.append(gwas_std_err)

	# Put everything in pandas dfs
	z_df = pd.DataFrame(np.transpose(z_score_mat), index=cafeh_studies, columns=cafeh_snps)
	n_df = pd.DataFrame(np.transpose(n_mat.astype(int)), index=cafeh_studies, columns=cafeh_snps)
	beta_df = pd.DataFrame(np.asmatrix(beta_arr), index=cafeh_studies[-1:], columns=cafeh_snps)
	std_err_df = pd.DataFrame(np.asmatrix(std_err_arr), index=cafeh_studies[-1:], columns=cafeh_snps)
	
	return z_df, beta_df, std_err_df, n_df

trait_name = sys.argv[1]  # Name of gwas trait
trait_sumstat_file = sys.argv[2]  # Summary statistic file for gwas trait
gtex_tissue_file = sys.argv[3]  # List of gtex tissues (and corresponding sample sizes)
genotype_reference_panel_dir = sys.argv[4]  # Directory containing genotype of variants
processed_gtex_associations_dir = sys.argv[5]  # Processed gtex associations (across all tissues)
processed_multivariate_cafeh_input_dir = sys.argv[6]  # output dir
chrom_num = int(sys.argv[7])
trait_sample_size = int(sys.argv[8])
cafeh_gene_list_file = sys.argv[9]



# Extract gtex tissue names and sample sizes from file
gtex_tissues, gtex_tissue_sample_sizes = get_gtex_tissue_names(gtex_tissue_file)


# Get an ordered list of genes on this chromosome as ell as a mapping to which tissues those genes are observed in
gene_to_tissues = get_genes_on_this_chromosome(cafeh_gene_list_file, 'chr' + str(chrom_num))
ordered_genes = np.asarray([*gene_to_tissues])


# Create a dictionary mapping from snp_id to gwas z_score (for snps on desired chromosome)
gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt = create_gwas_snp_id_to_beta_std_err_mapping_on_single_chromosome(trait_sumstat_file, str(chrom_num), trait_sample_size)


# Create output file handle
output_file = processed_multivariate_cafeh_input_dir + trait_name + '_chr' + str(chrom_num) + '_processed_gene_list.txt'
t = open(output_file,'w')
t.write('gene_id\tchrom_num\tgenotype_pkl_file\tz_pkl_file\tn_pkl_file\tbeta_pkl_file\tstd_err_pkl_file\n')



# Loop through genes
counter = 0
print(len(ordered_genes))
for gene_id in ordered_genes:
	counter = counter + 1
	print(counter)
	# Load in  processed_associations_dir data for this gene (contains info across all tissues)
	# Additionally filters out snps not present in all tissues
	gene_level_gtex_associations_file = processed_gtex_associations_dir + gene_id + '_hg19_associations.csv'
	gene_gtex_snps, gene_gtex_zscores, gene_gtex_tissues, gene_gtex_betas, gene_gtex_std_errs = extract_gene_level_gtex_data(gene_level_gtex_associations_file)

	# Get unique, ordered snps for this gene
	gene_gtex_unique_snps = np.sort(np.unique(gene_gtex_snps))

	# Get unique, ordered tissues for this gene
	cafeh_tissues = np.sort(np.unique(gene_gtex_tissues))

	# Filter gtex snps to those that also have gwas support
	cafeh_snps = filter_gtex_snps_to_those_that_have_gwas_support(gene_gtex_unique_snps, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt)

	# Skip genes with fewer than 10 snps
	if len(cafeh_snps) < 10:
		continue

	# HG38 to hg19 snp mapping file for this gene
	hg38_to_hg19_mapping_file = processed_gtex_associations_dir + gene_id + '_hg38_to_hg19_snp_mapping.txt'
	# Create mapping from hg38 snpid to hg19 snpid
	hg38_to_hg19_snp_id, hg19_to_hg38_snp_id = create_mapping_from_hg38_snp_id_to_hg19_snp_id(hg38_to_hg19_mapping_file)

	# Get genotype file (for computing r^2)
	genotype_file = genotype_reference_panel_dir + 'genotype_reference_panel_' + gene_id + '.txt'

	# convert genotype file to hg19, and order with snps
	genotype_df = convert_genotype_file_to_hg19_and_reorder(genotype_file, hg38_to_hg19_snp_id, hg19_to_hg38_snp_id, cafeh_snps)

	# Save genotype_df to pkl
	genotype_pkl_output_file = processed_multivariate_cafeh_input_dir + trait_name + '_' + gene_id + '_genotype_df.pkl'
	genotype_df.to_pickle(genotype_pkl_output_file)

	# Get summary stats in pandas data frame format
	z_df, beta_df, std_err_df, n_df = get_summary_stat_matrices(cafeh_snps, cafeh_tissues,  gene_gtex_snps, gene_gtex_tissues, gene_gtex_zscores, gene_gtex_betas, gene_gtex_std_errs, trait_name, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt)

	# Save summary stats to pickles
	# Save z_df to pkl
	z_pkl_output_file = processed_multivariate_cafeh_input_dir + trait_name + '_' + gene_id + '_z_df.pkl'
	z_df.to_pickle(z_pkl_output_file)
	# Save beta_df to pkl
	beta_pkl_output_file = processed_multivariate_cafeh_input_dir + trait_name + '_' + gene_id + '_beta_df.pkl'
	beta_df.to_pickle(beta_pkl_output_file)
	# Save std_err_df to pkl
	std_err_pkl_output_file = processed_multivariate_cafeh_input_dir + trait_name + '_' + gene_id + '_std_err_df.pkl'
	std_err_df.to_pickle(std_err_pkl_output_file)
	# Save n_df to pkl
	n_pkl_output_file = processed_multivariate_cafeh_input_dir + trait_name + '_' + gene_id + '_n_df.pkl'
	n_df.to_pickle(n_pkl_output_file)

	# Save locations of pickled file
	t.write(gene_id + '\t' + str(chrom_num) + '\t' + genotype_pkl_output_file + '\t' + z_pkl_output_file + '\t' + n_pkl_output_file + '\t' + beta_pkl_output_file + '\t' + std_err_pkl_output_file + '\n')


t.close()
