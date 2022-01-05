import numpy as np 
import os
import sys
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

def create_gwas_snp_id_to_zscore_mapping_on_single_chromosome(sumstat_file, chrom_string, trait_sample_size):
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
		# Extract relevent fields
		snp_id = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5]
		alt_snp_id = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4]

		# Get allele frequency
		af = float(data[6])
		info = float(data[7])
		beta_bolt_lmm_inf = float(data[10])
		se_bolt_lmm_inf = float(data[11])
		chi_sq_bolt_lmm = float(data[14])

		if chi_sq_bolt_lmm < 0:
			error_counter = error_counter +1
			continue
		unsigned_z_score = np.sqrt(chi_sq_bolt_lmm)

		if beta_bolt_lmm_inf < 0:
			zscore = unsigned_z_score*-1.0
		else:
			zscore = unsigned_z_score*1.0

		# Get maf from af
		if af > .5:
			maf = 1.0 - af
		else:
			maf = af
		if snp_id in dicti or alt_snp_id in alt_dicti:
			print('assumption oororor')
			pdb.set_trace()
		dicti[snp_id] = (zscore, maf, trait_sample_size, info)
		alt_dicti[alt_snp_id] = (zscore, maf, trait_sample_size, info)
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
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		# Quick error checking
		if len(data) != 11:
			continue
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[1]
		snp_id = data[3]
		maf = float(data[7])
		if maf > .5 or maf < .01:
			print('assumption eroror')
			pdb.set_trace()
		beta = float(data[9])
		std_err = float(data[10])
		zscore=beta/std_err
		snp_arr.append(snp_id)
		zscore_arr.append(zscore)
		tissue_arr.append(tissue_name)
	f.close()

	# Put in np arrays
	snp_arr = np.asarray(snp_arr)
	zscore_arr = np.asarray(zscore_arr)
	tissue_arr = np.asarray(tissue_arr)

	# Alphabetically order snps and zscores and tissues
	ordered_indices = np.argsort(snp_arr)

	return snp_arr[ordered_indices], zscore_arr[ordered_indices], tissue_arr[ordered_indices]

def get_bivariate_zscores_matrix_and_sample_size_matrix(tissue_snps, tissue_zscores, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt, gtex_tissue, trait_name, gtex_tissue_sample_size):
	bivariate_zscores = []
	snp_names = []
	bivariate_sample_sizes = []
	for snp_index, snp_name in enumerate(tissue_snps):
		tissue_zscore = tissue_zscores[snp_index]
		if snp_name in gwas_snp_id_to_z_score:
			if gwas_snp_id_to_z_score[snp_name][1] < .001 or gwas_snp_id_to_z_score[snp_name][3] < .6:
				continue
			bivariate_zscores.append((tissue_zscore, gwas_snp_id_to_z_score[snp_name][0]))
			snp_names.append(snp_name)
			bivariate_sample_sizes.append((gtex_tissue_sample_size, gwas_snp_id_to_z_score[snp_name][2]))
		elif snp_name in gwas_snp_id_to_z_score_alt:
			# ignore snp if maf < .1% or info score < .6
			if gwas_snp_id_to_z_score_alt[snp_name][1] < .001 or gwas_snp_id_to_z_score_alt[snp_name][3] < .6:
				continue
			bivariate_zscores.append((tissue_zscore, -1.0*gwas_snp_id_to_z_score_alt[snp_name][0]))
			snp_names.append(snp_name)
			bivariate_sample_sizes.append((gtex_tissue_sample_size, gwas_snp_id_to_z_score_alt[snp_name][2]))
		else:
			continue
	bivariate_zscore_mat = np.asarray(bivariate_zscores)
	snp_names = np.asarray(snp_names)
	bivariate_sample_size_mat = np.asarray(bivariate_sample_sizes)
	if len(snp_names) > 9:
		study_ids = np.asarray([gtex_tissue, trait_name])
		z_df = pd.DataFrame(np.transpose(bivariate_zscore_mat), index=study_ids, columns=snp_names)
		n_df = pd.DataFrame(np.transpose(bivariate_sample_size_mat.astype(int)), index=study_ids, columns=snp_names)
		boolers = True
	else:
		print('assumption error: Gene thrown out')
		z_df = 0
		n_df = 0
		snp_names = 0
		boolers = False
	return z_df, n_df, snp_names, boolers

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

def convert_genotype_file_to_hg19_and_reorder(karls_genotype_file, hg38_to_hg19_snp_id, hg19_to_hg38_snp_id, gene_gtex_unique_snps):
	karl_geno = np.loadtxt(karls_genotype_file, dtype=str)
	donor_id = karl_geno[1:,1]
	karl_geno = karl_geno[:,6:]
	hg38_snps = karl_geno[0,:]
	geno = karl_geno[1:,:]

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

	re_ordered_geno = mean_impute_missing_snps(re_ordered_geno)

	g_df = pd.DataFrame(re_ordered_geno, index=donor_id, columns=gene_gtex_unique_snps)

	return g_df

trait_name = sys.argv[1]  # Name of gwas trait
trait_sumstat_file = sys.argv[2]  # Summary statistic file for gwas trait
gtex_tissue_file = sys.argv[3]  # List of gtex tissues (and corresponding sample sizes)
gtex_cafeh_data = sys.argv[4]  # CAFEH Directory generated by karl
processed_gtex_associations_dir = sys.argv[5]  # Processed gtex associations (across all tissues)
processed_bivariate_cafeh_input_dir = sys.argv[6]  # output dir
chrom_num = int(sys.argv[7])
trait_sample_size = int(sys.argv[8])

# Cafeh gene list file
cafeh_gene_list_file = processed_gtex_associations_dir + 'cafeh_gene_list.txt'

# Extract gtex tissue names and sample sizes from file
gtex_tissues, gtex_tissue_sample_sizes = get_gtex_tissue_names(gtex_tissue_file)



# Get an ordered list of genes on this chromosome as ell as a mapping to which tissues those genes are observed in
gene_to_tissues = get_genes_on_this_chromosome(cafeh_gene_list_file, 'chr' + str(chrom_num))
ordered_genes = np.asarray([*gene_to_tissues])


# Create a dictionary mapping from snp_id to gwas z_score (for snps on desired chromosome)
gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt = create_gwas_snp_id_to_zscore_mapping_on_single_chromosome(trait_sumstat_file, str(chrom_num), trait_sample_size)


# Create array of output file handles
t_arr = []
for index, gtex_tissue in enumerate(gtex_tissues):
	output_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gtex_tissue + '_chr' + str(chrom_num) + '_processed_gene_list.txt'
	t_arr.append(open(output_file,'w'))
	t_arr[index].write('gene_id\tchrom_num\tgenotype_pkl_file\tz_pkl_file\tn_pkl_file\n')




# Loop through genes
counter = 0
print(len(ordered_genes))
for gene_id in ordered_genes:
	counter = counter + 1
	print(counter)
	# Load in  processed_associations_dir data for this gene (contains info across all tissues)
	gene_level_gtex_associations_file = processed_gtex_associations_dir + gene_id + '_hg19_associations.csv'
	gene_gtex_snps, gene_gtex_zscores, gene_gtex_tissues = extract_gene_level_gtex_data(gene_level_gtex_associations_file)
	# Get unique, ordered snps for this gene
	gene_gtex_unique_snps = np.sort(np.unique(gene_gtex_snps))
	# HG38 to hg19 snp mapping file for this gene
	hg38_to_hg19_mapping_file = processed_gtex_associations_dir + gene_id + '_hg38_to_hg19_snp_mapping.txt'
	# Create mapping from hg38 snpid to hg19 snpid
	hg38_to_hg19_snp_id, hg19_to_hg38_snp_id = create_mapping_from_hg38_snp_id_to_hg19_snp_id(hg38_to_hg19_mapping_file)

	# Get genotype file (for computing r^2)
	karls_genotype_file = gtex_cafeh_data + 'chr' + str(chrom_num) + '/' + gene_id + '/' + gene_id + '.raw'

	# convert genotype file to hg19, and order with snps
	genotype_df = convert_genotype_file_to_hg19_and_reorder(karls_genotype_file, hg38_to_hg19_snp_id, hg19_to_hg38_snp_id, gene_gtex_unique_snps)

	# Save genotype_df to pkl
	genotype_pkl_output_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gene_id + '_genotype_df.pkl'
	genotype_df.to_pickle(genotype_pkl_output_file)

	# loop through tisseus
	for tissue_index, gtex_tissue in enumerate(gtex_tissues):
		gtex_tissue_sample_size = gtex_tissue_sample_sizes[tissue_index]

		# Get indices corresponding to this tissue
		tissue_indices = gene_gtex_tissues == gtex_tissue
		# snps corresponding to this tissue
		tissue_snps = gene_gtex_snps[tissue_indices]
		tissue_zscores = gene_gtex_zscores[tissue_indices]

		if len(tissue_snps) < 10:
			gene_tissue_boolean = False
		else:
			# Get bivariate zscores for this tissues
			z_df, n_df, snp_names, gene_tissue_boolean = get_bivariate_zscores_matrix_and_sample_size_matrix(tissue_snps, tissue_zscores, gwas_snp_id_to_z_score, gwas_snp_id_to_z_score_alt, gtex_tissue, trait_name, gtex_tissue_sample_size)

		if gene_tissue_boolean == True:
			# Save z_df to pkl
			z_pkl_output_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gtex_tissue + '_' + gene_id + '_z_df.pkl'
			z_df.to_pickle(z_pkl_output_file)
			# Save n_df to pkl
			n_pkl_output_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gtex_tissue + '_' + gene_id + '_n_df.pkl'
			n_df.to_pickle(n_pkl_output_file)

			t_arr[tissue_index].write(gene_id + '\t' + str(chrom_num) + '\t' + genotype_pkl_output_file + '\t' + z_pkl_output_file + '\t' + n_pkl_output_file + '\n')

for index, gtex_tissue in enumerate(gtex_tissues):
	t_arr[index].close()
