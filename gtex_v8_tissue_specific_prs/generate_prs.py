import numpy as np 
import os
import sys
import pdb
from bgen.reader import BgenFile
import time

def create_mapping_from_full_samples_to_desired_samples(bgen_samples_file, desired_samples_file):
	full_samples = []
	desired_samples = []
	f = open(bgen_samples_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count < 2:
			head_count = head_count + 1
			continue
		if data[0] != data[1] or len(data) != 3:
			print('assumption eroror')
			pdb.set_trace()
		full_samples.append(data[0])
	f.close()

	f = open(desired_samples_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if len(data) != 2 or data[0] != data[1]:
			print('assumption eroror')
			pdb.set_trace()
		desired_samples.append(data[0])
	f.close()

	full_samples = np.asarray(full_samples)
	desired_samples = np.asarray(desired_samples)

	full_sample_to_position = {}
	for i, full_sample in enumerate(full_samples):
		full_sample_to_position[full_sample] = i

	mapping = []

	for desired_sample in desired_samples:
		mapping.append(full_sample_to_position[desired_sample])

	mapping = np.asarray(mapping)

	if np.array_equal(full_samples[mapping], desired_samples) == False:
		print('assumpton eororor')

	return mapping, len(full_samples), len(desired_samples), desired_samples

def get_tissue_names_from_prs_betas_file(cafeh_prs_betas_file):
	f = open(cafeh_prs_betas_file)
	head_count = 0
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			tissue_names = np.asarray(data[1:])
			continue
		break
	f.close()
	return tissue_names

def create_mapping_from_variant_to_effect_sizes(cafeh_prs_betas_file):
	f = open(cafeh_prs_betas_file)
	head_count = 0
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			tissue_names = np.asarray(data[1:])
			continue
		variant_id = data[0]
		weights = np.asarray(data[1:]).astype(float)
		# Fillter betas
		#for index, weight in enumerate(weights):
			#if np.abs(weight) > beta_threshold:
				#weights[index] = 0.0
		# Done fillter betas
		variant_info = data[0].split('_')
		chrom_num = variant_info[0].split('hr')[1]

		reformated_variant_id = chrom_num + ':' + variant_info[1] + '_' + variant_info[2] + '_' + variant_info[3]
		reformated_variant_id_alt = chrom_num + ':' + variant_info[1] + '_' + variant_info[3] + '_' + variant_info[2]

		if reformated_variant_id in mapping or reformated_variant_id_alt in mapping:
			print('assumption eroor')
			continue
		mapping[reformated_variant_id] = np.reshape(weights, (1,len(weights)))
		mapping[reformated_variant_id_alt] = -1.0*np.reshape(weights, (1,len(weights)))
	f.close()
	return mapping


bgen_file = sys.argv[1]
bgen_samples_file = sys.argv[2]
desired_samples_file = sys.argv[3]
cafeh_prs_betas_dir = sys.argv[4] # input file
ukbb_prs_dir = sys.argv[5] # input file
chrom = sys.argv[6]
trait_name = sys.argv[7]



# We are going create seperate prs's iterating over the following parameters
coloc_thresholds = [.1, .3, .5, .7, .9]

methods = ['coloc', 'causal_v1_coloc', 'causal_v2_coloc', 'mm_v1_coloc', 'mm_v2_coloc']



versions = []
variant_to_effect_sizes_arr = []
sizes = []
all_variants = {}
for coloc_threshold in coloc_thresholds:
	for method in methods:

		prs_output_file = ukbb_prs_dir + trait_name + '_' + str(coloc_threshold) + '_' + method + '_prs_chrom_' + chrom + '.txt'
		#prs_beta_file = cafeh_prs_betas_dir + 'cafeh_results_' + trait_name + '_' + method + '_' + str(p_active_threshold) + '_prs_beta_chrom_' + chrom + '.txt'
		prs_beta_file = cafeh_prs_betas_dir + method + '_results_' + trait_name + '_' + str(coloc_threshold) + '_prs_beta_chrom_' + chrom + '.txt'
		versions.append((str(coloc_threshold) + '_' + method, prs_output_file))
		variant_to_effect_sizes = create_mapping_from_variant_to_effect_sizes(prs_beta_file)
		print(len(variant_to_effect_sizes))
		sizes.append(len(variant_to_effect_sizes))
		variant_to_effect_sizes_arr.append(variant_to_effect_sizes)
		variants = [*variant_to_effect_sizes]
		for variant in variants:
			all_variants[variant] = 1


tissue_names = get_tissue_names_from_prs_betas_file(prs_beta_file)


sample_mapping, total_samples, num_prs_samples, prs_sample_names = create_mapping_from_full_samples_to_desired_samples(bgen_samples_file, desired_samples_file)


# Create matrix to keep track of PRS
prs_mats = []
for version in versions:
	prs_mat = np.zeros((num_prs_samples, len(tissue_names)))
	prs_mats.append(prs_mat)


counter = 1

used = {}

with BgenFile(bgen_file, delay_parsing=True) as bfile:
	for var in bfile:
		dosage = var.minor_allele_dosage
		temp_variant_id = var.varid
		minor_allele = var.minor_allele
		alleles = var.alleles
		# error checking
		if len(dosage) != total_samples:
			print('assumption eroror')
		if len(alleles) != 2:
			print('assumption eororro')

		allele = 'False'
		for allele in alleles:
			if allele != minor_allele:
				major_allele = allele
		if allele == 'False':
			print('assumption erroror')
			pdb.set_trace()
		
		# Want variant id to be of form chrom_num:position_non-effect-allele_effect_allele
		variant_id = temp_variant_id.split('_')[0] + '_' + major_allele + '_' + minor_allele

		variant_info = variant_id.split('_')
		if len(variant_info) != 3:
			variant_id = var.chrom + ':' + str(var.pos) + '_' + var.alleles[0] + '_' + var.alleles[1]
			variant_info = variant_id.split('_')
		counter = counter + 1

		if variant_info[2] != minor_allele:
			print('minor allele assumption erororor')
			pdb.set_trace()
		if variant_id not in all_variants:
			continue
		########################################
		# SHOULD ONLY DO THIS STUFF IF PASS DICTI
		###########################################
		# convert dosage to dosage for desired samples only
		dosage_desired_samples = dosage[sample_mapping]

		dosage_desired_samples = np.reshape(dosage_desired_samples,(len(dosage_desired_samples),1))

		#mid_time = time.time()
		used[variant_id] = 1
		# Get prs effect sizes
		for version_index, version in enumerate(versions):
			if variant_id not in variant_to_effect_sizes_arr[version_index]:
				continue
			effect_sizes = variant_to_effect_sizes_arr[version_index][variant_id]
			prs_mats[version_index] = prs_mats[version_index] + np.dot(dosage_desired_samples, effect_sizes)
			#for tissue_index in range(len(tissue_names)):
			#	prs_mats[version_index][:, tissue_index] = prs_mats[version_index][:, tissue_index] + (dosage_desired_samples*effect_sizes[tissue_index])


if len(all_variants) > 0:
	variant_efficiency = len(used)/(len(all_variants)/2.0)
	print('variant efficiency: ' + str(variant_efficiency))


header = []
header.append('sample_name')
for tissue in tissue_names:
	header.append(tissue)

for version_index, version_tuple in enumerate(versions):

	temp = np.hstack((np.asmatrix(prs_sample_names).T, prs_mats[version_index].astype(str)))

	final_mat = np.vstack((np.asmatrix(header), temp))

	np.savetxt(version_tuple[1], final_mat, fmt="%s", delimiter='\t')

