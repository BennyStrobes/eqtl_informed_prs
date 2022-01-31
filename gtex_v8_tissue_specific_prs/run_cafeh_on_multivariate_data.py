import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import cafeh
import pandas as pd
from cafeh.cafeh_summary import fit_cafeh_summary, fit_cafeh_z, fit_cafeh_summary_fixed_pi, fit_cafeh_z_fixed_pi
from cafeh.model_queries import *


def filter_snps_in_g_df(df1, df2):
	return df1[df2.columns.values]


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
	
	return colocalized_variant_report, colocalized_components_dicti

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


def get_cafeh_predicted_snp_eqtl_p_active(snp_names, colocalized_variant_report, tissue_name):
	trait_colocalized_variant_report = colocalized_variant_report[colocalized_variant_report['study'] == tissue_name]

	snp_predicted_effects_dicti = {}
	coloc_per_component_effect_size = np.asarray(trait_colocalized_variant_report['p_active'])
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



def print_snp_predicted_effects_for_a_gene(snp_names, ordered_tissues, pred_effects, snp_predicted_effects_file):

	snp_predicted_effects_str = pred_effects.astype(str)

	t = open(snp_predicted_effects_file,'w')
	t.write('snp_id\t' + '\t'.join(ordered_tissues) + '\n')

	for snp_index, snp_name in enumerate(snp_names):
		snp_predicted_effect = snp_predicted_effects_str[snp_index,:]
		t.write(snp_name + '\t' + '\t'.join(snp_predicted_effect) + '\n')
	t.close()

def print_num_eqtl_components(num_eqtl_components, num_eqtl_components_output_file):
	t = open(num_eqtl_components_output_file,'w')
	t.write(str(num_eqtl_components) + '\n')
	t.close()

def extract_component_predicted_effects_from_cafeh_model(cafeh_fixed_model, colocalizing_component):
	# fixed_variant_report = summary_table(cafeh_fixed_model, filter_variants=False, min_p_active=0.0)
	# fixed_variant_report[fixed_variant_report['top_component'] == colocalizing_component]
	tmp_component_effect = cafeh_fixed_model.weight_means[0,colocalizing_component,:]
	tmp_p_active = cafeh_fixed_model.active[0,colocalizing_component]
	tmp_pi = cafeh_fixed_model.pi[colocalizing_component,:]
	# Prune pi 
	tmp_pi[tmp_pi < .001] = 0.0

	# Get predicted effect
	predicted_effect = tmp_pi*tmp_p_active*tmp_component_effect
	return predicted_effect

def print_num_eqtl_components(cafeh_z_model, ordered_tissues, num_eqtl_components_output_file):
	cafeh_studies = cafeh_z_model.study_ids
	components_per_study = np.sum(cafeh_z_model.active > .5,axis=1)
	mapping = {}
	for index, cafeh_study in enumerate(cafeh_studies):
		mapping[cafeh_study] = components_per_study[index]

	t = open(num_eqtl_components_output_file,'w')
	t.write('tissue\tnum_components\n')
	for tissue in ordered_tissues:
		if tissue not in mapping:
			counter = 0
		else:
			counter = mapping[tissue]
		t.write(tissue + '\t' + str(counter) + '\n')
	t.close()

def print_num_coloc_components(colocalizing_components, ordered_tissues, num_coloc_components_output_file):
	mapping = {}
	for tissue in ordered_tissues:
		mapping[tissue] = 0

	for colocalizing_component in colocalizing_components:
		colocalizing_tissues = colocalizing_component[2]
		for colocalizing_tissue in colocalizing_tissues:
			mapping[colocalizing_tissue] = mapping[colocalizing_tissue] + 1

	t = open(num_coloc_components_output_file,'w')
	t.write('tissue\tnum_components\n')
	for tissue in ordered_tissues:
		counter = mapping[tissue]
		t.write(tissue + '\t' + str(counter) + '\n')
	t.close()	

def cafeh_wrapper_debug(gene_name, genotype_file, zscore_file, sample_size_file, beta_file, std_err_file, trait_name, ordered_tissues, tissue_to_position_mapping, output_stem, version):
	# Load in pickled cafeh input data
	g_df = pd.read_pickle(genotype_file)
	z_df = pd.read_pickle(zscore_file)
	n_df = pd.read_pickle(sample_size_file)
	beta_df = pd.read_pickle(beta_file)
	std_err_df = pd.read_pickle(std_err_file)


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
	#g_df = g_df.astype(float)
	LD = np.corrcoef(np.transpose(g_df))
	LD_df = pd.DataFrame(LD, index=g_df.columns.values, columns=g_df.columns.values)

	#######################
	# Run CAFEH
	#######################
	cafeh_z_model = fit_cafeh_z(LD_df, z_df, n=n_df,K=10, prior_variance=10.0)


	# get names of active cafeh trait components components
	active_trait_components = np.where(cafeh_z_model.active[-1,:] > .95)[0]

	if len(active_trait_components) > 5:
		active_trait_components = np.asarray([])

	# Get names of colocalizing (with any tissue) cafeh components
	colocalizing_components = []
	for active_trait_component in active_trait_components:
		if np.sum(cafeh_z_model.active[:,active_trait_component] > .95) > 1:
			tissue_indices = np.where(cafeh_z_model.active[:-1,active_trait_component] > .95)[0]
			p_active_arr = cafeh_z_model.active[tissue_indices,active_trait_component]
			tissue_name_arr = cafeh_z_model.study_ids[tissue_indices]
			colocalizing_components.append((active_trait_component, p_active_arr, tissue_name_arr))
			print(tissue_name_arr)

	'''
	# Initialize matrix of predicted effects
	pred_effects = np.zeros((len(snp_names), len(ordered_tissues)))
	weighted_pred_effects = np.zeros((len(snp_names), len(ordered_tissues)))

	# IF there are colocalizing snps, run cafeh again with pi's fixed using beta, std-err model
	if len(colocalizing_components) > 0:
		subset_n_df = n_df.iloc[-1:,:]
		cafeh_fixed_model = fit_cafeh_summary_fixed_pi(LD_df, beta_df, std_err_df, np.copy(cafeh_z_model.pi), n=subset_n_df, K=10)

		# Loop through colocalizing components
		for colocalizing_component_tuple in colocalizing_components:
			# Get predicted trait effect from this component
			colocalizing_component = colocalizing_component_tuple[0]
			component_predicted_effects = extract_component_predicted_effects_from_cafeh_model(cafeh_fixed_model, colocalizing_component)

			# Get tissues colocalizing with this component
			colocalizing_tissues = colocalizing_component_tuple[2]
			colocalizing_p_actives = colocalizing_component_tuple[1]

			# Loop through colocalizing tissues
			for tissue_index, colocalizing_tissue in enumerate(colocalizing_tissues):
				colocalizing_p_active = colocalizing_p_actives[tissue_index]
				column_index = tissue_to_position_mapping[colocalizing_tissue]

				pred_effects[:, column_index] = pred_effects[:, column_index] + component_predicted_effects
				weighted_pred_effects[:, column_index] = weighted_pred_effects[:, column_index] + (colocalizing_p_active*component_predicted_effects)
		
		# Save some files if we get colocalization for this gene
		p_active_file = output_stem + gene_name + '_p_active.txt'
		np.savetxt(p_active_file, np.hstack((np.asmatrix(cafeh_z_model.study_ids).T, cafeh_z_model.active.astype(str))), fmt="%s", delimiter='\t')
		pi_file = output_stem + gene_name + '_pi.txt'
		np.savetxt(pi_file, np.vstack((np.asmatrix(cafeh_z_model.snp_ids), cafeh_z_model.pi.astype(str))), fmt="%s", delimiter='\t')

	
	# print predicted effects of snps for this gene
	snp_predicted_effects_file = output_stem + gene_name + '_predicted_effects.txt'
	print_snp_predicted_effects_for_a_gene(snp_names, ordered_tissues, pred_effects, snp_predicted_effects_file)
	weighted_snp_predicted_effects_file = output_stem + gene_name + '_weighted_predicted_effects.txt'
	print_snp_predicted_effects_for_a_gene(snp_names, ordered_tissues, weighted_pred_effects, weighted_snp_predicted_effects_file)
	
	# print number of eqtl components for this gene
	num_eqtl_components_output_file = output_stem + gene_name + '_number_eqtl_components.txt'
	print_num_eqtl_components(cafeh_z_model, ordered_tissues, num_eqtl_components_output_file)

	# print number of colocalizing components
	num_coloc_components_output_file = output_stem + gene_name + '_number_coloc_components.txt'
	print_num_coloc_components(colocalizing_components, ordered_tissues, num_coloc_components_output_file)
	'''

def cafeh_wrapper(gene_name, genotype_file, zscore_file, sample_size_file, beta_file, std_err_file, trait_name, ordered_tissues, tissue_to_position_mapping, output_stem, version):
	# Load in pickled cafeh input data
	g_df = pd.read_pickle(genotype_file)
	z_df = pd.read_pickle(zscore_file)
	n_df = pd.read_pickle(sample_size_file)
	beta_df = pd.read_pickle(beta_file)
	std_err_df = pd.read_pickle(std_err_file)


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
	#g_df = g_df.astype(float)
	LD = np.corrcoef(np.transpose(g_df))
	LD_df = pd.DataFrame(LD, index=g_df.columns.values, columns=g_df.columns.values)

	#######################
	# Run CAFEH
	#######################
	cafeh_z_model = fit_cafeh_z(LD_df, z_df, n=n_df,K=10, prior_variance=10.0)


	if version.startswith('all_tissues_.'):
		p_active_thresh = float(version.split('_')[2])
	else:
		p_active_thresh = .5
	# get names of active cafeh trait components components
	active_trait_components = np.where(cafeh_z_model.active[-1,:] > p_active_thresh)[0]


	# Get names of colocalizing (with any tissue) cafeh components
	colocalizing_components = []
	for active_trait_component in active_trait_components:
		if np.sum(cafeh_z_model.active[:,active_trait_component] > p_active_thresh) > 1:
			if version.startswith('all_tissues'):
				tissue_indices = np.where(cafeh_z_model.active[:-1,active_trait_component] > p_active_thresh)[0]
			elif version == 'top_tissue':
				tissue_indices = np.asarray([np.argmax(cafeh_z_model.active[:-1,active_trait_component])])
			else:
				print('assumptino erroror')
				pdb.set_trace()
			p_active_arr = cafeh_z_model.active[tissue_indices,active_trait_component]
			tissue_name_arr = cafeh_z_model.study_ids[tissue_indices]
			print(tissue_name_arr)
			colocalizing_components.append((active_trait_component, p_active_arr, tissue_name_arr))


	# Initialize matrix of predicted effects
	pred_effects = np.zeros((len(snp_names), len(ordered_tissues)))
	weighted_pred_effects = np.zeros((len(snp_names), len(ordered_tissues)))

	# IF there are colocalizing snps, run cafeh again with pi's fixed using beta, std-err model
	if len(colocalizing_components) > 0:
		subset_n_df = n_df.iloc[-1:,:]
		cafeh_fixed_model = fit_cafeh_summary_fixed_pi(LD_df, beta_df, std_err_df, np.copy(cafeh_z_model.pi), n=subset_n_df, K=10)

		# Loop through colocalizing components
		for colocalizing_component_tuple in colocalizing_components:
			# Get predicted trait effect from this component
			colocalizing_component = colocalizing_component_tuple[0]
			component_predicted_effects = extract_component_predicted_effects_from_cafeh_model(cafeh_fixed_model, colocalizing_component)

			# Get tissues colocalizing with this component
			colocalizing_tissues = colocalizing_component_tuple[2]
			colocalizing_p_actives = colocalizing_component_tuple[1]

			# Loop through colocalizing tissues
			for tissue_index, colocalizing_tissue in enumerate(colocalizing_tissues):
				colocalizing_p_active = colocalizing_p_actives[tissue_index]
				column_index = tissue_to_position_mapping[colocalizing_tissue]

				pred_effects[:, column_index] = pred_effects[:, column_index] + component_predicted_effects
				weighted_pred_effects[:, column_index] = weighted_pred_effects[:, column_index] + (colocalizing_p_active*component_predicted_effects)
		
		# Save some files if we get colocalization for this gene
		p_active_file = output_stem + gene_name + '_p_active.txt'
		np.savetxt(p_active_file, np.hstack((np.asmatrix(cafeh_z_model.study_ids).T, cafeh_z_model.active.astype(str))), fmt="%s", delimiter='\t')
		pi_file = output_stem + gene_name + '_pi.txt'
		np.savetxt(pi_file, np.vstack((np.asmatrix(cafeh_z_model.snp_ids), cafeh_z_model.pi.astype(str))), fmt="%s", delimiter='\t')

	
	# print predicted effects of snps for this gene
	snp_predicted_effects_file = output_stem + gene_name + '_predicted_effects.txt'
	print_snp_predicted_effects_for_a_gene(snp_names, ordered_tissues, pred_effects, snp_predicted_effects_file)
	weighted_snp_predicted_effects_file = output_stem + gene_name + '_weighted_predicted_effects.txt'
	print_snp_predicted_effects_for_a_gene(snp_names, ordered_tissues, weighted_pred_effects, weighted_snp_predicted_effects_file)
	
	# print number of eqtl components for this gene
	num_eqtl_components_output_file = output_stem + gene_name + '_number_eqtl_components.txt'
	print_num_eqtl_components(cafeh_z_model, ordered_tissues, num_eqtl_components_output_file)

	# print number of colocalizing components
	num_coloc_components_output_file = output_stem + gene_name + '_number_coloc_components.txt'
	print_num_coloc_components(colocalizing_components, ordered_tissues, num_coloc_components_output_file)



def get_num_genes(gene_file):
	f = open(gene_file)
	gene_counter = 0
	for line in f:
		gene_counter = gene_counter + 1
	f.close()
	return (gene_counter-1)

# For parallelization purposes
def parallelization_start_and_end(num_tasks, job_number, total_jobs):
	tasks_per_job = (num_tasks/total_jobs) + 1
	start_task = job_number*tasks_per_job
	end_task = (job_number + 1)*tasks_per_job -1 
	return start_task, end_task

def get_ordered_tissues(gtex_tissue_file):
	f = open(gtex_tissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	arr = np.asarray(arr)
	dicti = {}
	for i, ele in enumerate(arr):
		dicti[ele] = i
	return arr, dicti


gene_file = sys.argv[1]
trait_name = sys.argv[2]
gtex_tissue_file = sys.argv[3]
cafeh_output_dir = sys.argv[4]
job_number = int(sys.argv[5])
total_jobs = int(sys.argv[6])
version = sys.argv[7]


print(trait_name)
print(job_number)

# Extract names of gtex tissues
ordered_tissues, tissue_to_position_mapping = get_ordered_tissues(gtex_tissue_file)



#For parallelization purposes
num_genes = get_num_genes(gene_file)
start_number, end_number = parallelization_start_and_end(num_genes, job_number, total_jobs)

# file stem to write results to
output_stem = cafeh_output_dir + 'cafeh_results_' + trait_name + '_' + version + '_'

head_count = 0
counter = -2
f = open(gene_file)
for line in f:
	counter = counter + 1
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Skip gene ids not in this parallelization run
	if counter < start_number or counter > end_number:
		continue

	# Extract relevent fields
	gene_name = data[0]
	print(gene_name)

	gene_chrom = data[1]
	genotype_file = data[2]
	zscore_file = data[3]
	sample_size_file = data[4]
	beta_file = data[5]
	std_err_file = data[6]
	if version != 'all_tissues_debug':
		cafeh_wrapper(gene_name, genotype_file, zscore_file, sample_size_file, beta_file, std_err_file, trait_name, ordered_tissues, tissue_to_position_mapping, output_stem, version)
	else:
		cafeh_wrapper_debug(gene_name, genotype_file, zscore_file, sample_size_file, beta_file, std_err_file, trait_name, ordered_tissues, tissue_to_position_mapping, output_stem, version)
f.close()
