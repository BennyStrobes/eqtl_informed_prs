import numpy as np 
import os
import sys
import pdb
import cafeh
import pandas as pd
from cafeh.cafeh_summary import fit_cafeh_summary, fit_cafeh_z
from cafeh.model_queries import *


def filter_snps_in_g_df(df1, df2):
	return df1[df2.columns.values]

def debugging(g_df, gene_name):
	snp_hg38_to_hg19_mapping_file = '/work-zfs/abattle4/bstrober/eqtl_informed_prs/gtex_v8/processed_gtex_associations/' + gene_name + '_hg38_to_hg19_snp_mapping.txt'
	hg19_to_hg38 = {}
	hg38_to_hg19 = {}
	f = open(snp_hg38_to_hg19_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		hg38 = data[0]
		hg19 = data[1]
		if hg19 == 'NA':
			continue
		hg19_to_hg38[hg19] = hg38
		hg38_to_hg19[hg38] = hg19
	f.close()

	chrom_string = g_df.columns[0].split('_')[0]

	raw_hg38_genotype_file = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/' + chrom_string + '/' + gene_name + '/' + gene_name + '.raw'

	raw_hg38_genotype = np.loadtxt(raw_hg38_genotype_file,dtype=str,delimiter=' ')
	raw_hg38_genotype = raw_hg38_genotype[:,6:]

	num_snps = raw_hg38_genotype.shape[1]

	for snp_num in range(num_snps):
		hg38_snp_name_raw = raw_hg38_genotype[0, snp_num]
		hg38_snp_name = hg38_snp_name_raw.split('_b38')[0] + '_b38'
		info = hg38_snp_name_raw.split('_')
		if len(info) != 6:
			print('assuption error')
			pdb.set_trace()
		if info[3] != info[5]:
			print('assusmptionerororo')
			pdb.set_trace()
		if hg38_snp_name not in hg38_to_hg19:
			continue
		genotype = raw_hg38_genotype[1:, snp_num]
		observed_indices = genotype != 'NA'
		genotype = genotype[observed_indices]
		hg19_snp_name = hg38_to_hg19[hg38_snp_name]
		if hg19_snp_name not in g_df.columns:
			continue
		corry = np.corrcoef(np.asarray(g_df[hg19_snp_name])[observed_indices], genotype.astype(float))
		if corry[0,1] < .99:
			print('erroror')
			pdb.set_trace()
		#if corry[0,1] != 1.0:
		if np.array_equal(np.asarray(g_df[hg19_snp_name])[observed_indices], genotype.astype(float)) == False:
			print('eroror')
			pdb.set_trace()

def get_colocalized_summary_table(variant_report, trait_name, tissue_name):
	# Get components active in gwas study
	trait_variant_report = variant_report[variant_report['study'] == trait_name]
	trait_components = np.unique(np.asarray(trait_variant_report[trait_variant_report['p_active'] > .5]['top_component']))

	# Get components active in eqtl study
	tissue_variant_report = variant_report[variant_report['study'] == tissue_name]
	tissue_components = np.unique(np.asarray(tissue_variant_report[tissue_variant_report['p_active'] > .5]['top_component']))

	# Get colocalizing components
	colocalized_components_dicti = {}
	for trait_component in trait_components:
		if trait_component in tissue_components:
			colocalized_components_dicti[trait_component] = 1
	
	# Extract colocalized variant report
	colocalized_variant_report = variant_report[variant_report['top_component'].isin(colocalized_components_dicti)]
	
	return colocalized_variant_report

def get_cafeh_predicted_snp_effects_from_colocalizing_components(snp_names, colocalized_variant_report, trait_name):
	trait_colocalized_variant_report = colocalized_variant_report[colocalized_variant_report['study'] == trait_name]
	
	snp_predicted_effects_dicti = {}
	coloc_per_component_effect_size = np.asarray(trait_colocalized_variant_report['p_active']*trait_colocalized_variant_report['effect']*trait_colocalized_variant_report['pi'])
	coloc_ordered_variants = np.asarray(trait_colocalized_variant_report['variant_id'])
	for indexer, variant_id in enumerate(coloc_ordered_variants):
		snp_predicted_effects_dicti[variant_id] = coloc_per_component_effect_size[indexer]

	snp_predicted_effects = []
	for snp_name in snp_names:
		if snp_name in snp_predicted_effects_dicti:
			snp_predicted_effects.append(snp_predicted_effects_dicti[snp_name])
		else:
			snp_predicted_effects.append(0.0)
	return np.asarray(snp_predicted_effects)

def cafeh_wrapper(gene_name, genotype_file, zscore_file, sample_size_file, trait_name, tissue_name):
	# Load in pickled cafeh input data
	g_df = pd.read_pickle(genotype_file)
	z_df = pd.read_pickle(zscore_file)
	n_df = pd.read_pickle(sample_size_file)

	# Need to filter g_df to only contain columns found in z
	g_df = filter_snps_in_g_df(g_df, z_df)

	# Extract snp names
	snp_names = np.asarray(g_df.columns)

	# Quick error checking
	if np.array_equal(g_df.columns.values, z_df.columns.values) == False:
		print('assumption oerooror')
		pdb.set_trace()
	if np.array_equal(z_df.columns.values, n_df.columns.values) == False:
		print('assumption eororor')
		pdb.set_trace()

	# Convert from genotype space to LD space
	LD = np.corrcoef(np.transpose(g_df))
	LD_df = pd.DataFrame(LD, index=g_df.columns.values, columns=g_df.columns.values)

	#######################
	# Run CAFEH
	#######################
	cafehs = fit_cafeh_z(LD_df, z_df, n=n_df,K=10, prior_variance=10.0)
	# Get un-filtered variant report
	variant_report = summary_table(cafehs, filter_variants=False, min_p_active=0.0, max_snps= (LD.shape[0])*2 + 1)
	# Filter to colocalized variants
	colocalized_variant_report = get_colocalized_summary_table(variant_report, trait_name, tissue_name)
	# Get predicted effect from all snps (limiting to colocalizing components)
	snp_predicted_effects = get_cafeh_predicted_snp_effects_from_colocalizing_components(snp_names, colocalized_variant_report, trait_name)

	#######################
	# print to output file
	#######################	

gene_file = sys.argv[1]
trait_name = sys.argv[2]
tissue_name = sys.argv[3]
bivariate_cafeh_output_dir = sys.argv[4]



head_count = 0
counter = 0
print(gene_file)
f = open(gene_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	gene_name = data[0]
	print(gene_name)

	gene_chrom = data[1]
	genotype_file = data[2]
	zscore_file = data[3]
	sample_size_file = data[4]
	cafeh_wrapper(gene_name, genotype_file, zscore_file, sample_size_file, trait_name, tissue_name)


	counter = counter + 1