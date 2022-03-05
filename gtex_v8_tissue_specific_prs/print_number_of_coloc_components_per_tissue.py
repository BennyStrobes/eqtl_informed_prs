import numpy as np 
import os
import sys
import pdb



def get_avg_number_of_expressed_genes_across_composit_tissues(composit_tissues):
	num_genes = []
	for composit_tissue in composit_tissues:
		expression_file = '/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/' + composit_tissue + '_normalized_expression_autosomes_protein_coding_and_linc.txt'
		tissue_count = -1
		f = open(expression_file)
		for line in f:
			tissue_count = tissue_count + 1
		f.close()
		num_genes.append(tissue_count)
	num_genes = np.asarray(num_genes)
	return np.mean(num_genes)

gtex_tissue_file = sys.argv[1]
trait_name = sys.argv[2]
bivariate_cafeh_output_dir = sys.argv[3]



tissues = []
tissue_to_sample_size = {}
tissue_to_num_components = {}
tissue_to_num_eqtls = {}
tissue_to_number_expressed_genes = {}


f = open(gtex_tissue_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	samp_size = data[1]
	tissues.append(tissue_name)
	tissue_to_sample_size[tissue_name] = samp_size
	tissue_to_num_components[tissue_name] = 0
	tissue_to_num_eqtls[tissue_name] = 0

	composit_tissues = data[3].split(',')
	avg_num_genes = get_avg_number_of_expressed_genes_across_composit_tissues(composit_tissues)
	tissue_to_number_expressed_genes[tissue_name] = avg_num_genes

f.close()
tissues = np.asarray(tissues)


coloc_thresholds = [.1, .3, .5, .7, .9]

for coloc_threshold in coloc_thresholds:
	for tissue in tissues:
		tissue_to_num_components[tissue] = 0
		tissue_to_num_eqtls[tissue] = 0

	for chrom_num in range(1,23):
		component_file = bivariate_cafeh_output_dir + 'coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

		f = open(component_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			num_components = int(data[1])
			num_eqtls = int(data[2])
			tissue_to_num_components[tissue_name] = tissue_to_num_components[tissue_name] + num_components
			tissue_to_num_eqtls[tissue_name] = tissue_to_num_eqtls[tissue_name] + num_eqtls
		f.close()

	# print to output file
	t = open(bivariate_cafeh_output_dir + 'coloc_results_' + trait_name+ '_' + str(coloc_threshold) + '_num_prs_components.txt', 'w')
	t.write('tissue\tsample_size\tnum_cafeh_components\tnum_eqtl_components\tavg_num_expressed_genes\n')
	for tissue in tissues:
		samp_size = tissue_to_sample_size[tissue]
		num_comp = tissue_to_num_components[tissue]
		num_eqtls = tissue_to_num_eqtls[tissue]
		avg_num_expressed_genes = tissue_to_number_expressed_genes[tissue]
		t.write(tissue + '\t' + str(samp_size) + '\t' + str(num_comp) + '\t' + str(num_eqtls) + '\t' + str(avg_num_expressed_genes) + '\n')
	t.close()










tissues = []
tissue_to_sample_size = {}
tissue_to_num_components = {}
tissue_to_num_eqtls = {}
tissue_to_number_expressed_genes = {}


f = open(gtex_tissue_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	samp_size = data[1]
	tissues.append(tissue_name)
	tissue_to_sample_size[tissue_name] = samp_size
	tissue_to_num_components[tissue_name] = 0
	tissue_to_num_eqtls[tissue_name] = 0

	composit_tissues = data[3].split(',')
	avg_num_genes = get_avg_number_of_expressed_genes_across_composit_tissues(composit_tissues)
	tissue_to_number_expressed_genes[tissue_name] = avg_num_genes

f.close()
tissues = np.asarray(tissues)


coloc_thresholds = [.1, .3, .5, .7, .9]

for coloc_threshold in coloc_thresholds:
	for tissue in tissues:
		tissue_to_num_components[tissue] = 0
		tissue_to_num_eqtls[tissue] = 0

	for chrom_num in range(1,23):
		component_file = bivariate_cafeh_output_dir + 'causal_v1_coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

		f = open(component_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			num_components = int(data[1])
			num_eqtls = int(data[2])
			tissue_to_num_components[tissue_name] = tissue_to_num_components[tissue_name] + num_components
			tissue_to_num_eqtls[tissue_name] = tissue_to_num_eqtls[tissue_name] + num_eqtls
		f.close()

	# print to output file
	t = open(bivariate_cafeh_output_dir + 'causal_v1_coloc_results_' + trait_name+ '_' + str(coloc_threshold) + '_num_prs_components.txt', 'w')
	t.write('tissue\tsample_size\tnum_cafeh_components\tnum_eqtl_components\tavg_num_expressed_genes\n')
	for tissue in tissues:
		samp_size = tissue_to_sample_size[tissue]
		num_comp = tissue_to_num_components[tissue]
		num_eqtls = tissue_to_num_eqtls[tissue]
		avg_num_expressed_genes = tissue_to_number_expressed_genes[tissue]
		t.write(tissue + '\t' + str(samp_size) + '\t' + str(num_comp) + '\t' + str(num_eqtls) + '\t' + str(avg_num_expressed_genes) + '\n')
	t.close()









tissues = []
tissue_to_sample_size = {}
tissue_to_num_components = {}
tissue_to_num_eqtls = {}
tissue_to_number_expressed_genes = {}


f = open(gtex_tissue_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	samp_size = data[1]
	tissues.append(tissue_name)
	tissue_to_sample_size[tissue_name] = samp_size
	tissue_to_num_components[tissue_name] = 0
	tissue_to_num_eqtls[tissue_name] = 0

	composit_tissues = data[3].split(',')
	avg_num_genes = get_avg_number_of_expressed_genes_across_composit_tissues(composit_tissues)
	tissue_to_number_expressed_genes[tissue_name] = avg_num_genes

f.close()
tissues = np.asarray(tissues)


coloc_thresholds = [.1, .3, .5, .7, .9]

for coloc_threshold in coloc_thresholds:
	for tissue in tissues:
		tissue_to_num_components[tissue] = 0
		tissue_to_num_eqtls[tissue] = 0

	for chrom_num in range(1,23):
		component_file = bivariate_cafeh_output_dir + 'causal_v2_coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

		f = open(component_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			num_components = int(data[1])
			num_eqtls = int(data[2])
			tissue_to_num_components[tissue_name] = tissue_to_num_components[tissue_name] + num_components
			tissue_to_num_eqtls[tissue_name] = tissue_to_num_eqtls[tissue_name] + num_eqtls
		f.close()

	# print to output file
	t = open(bivariate_cafeh_output_dir + 'causal_v2_coloc_results_' + trait_name+ '_' + str(coloc_threshold) + '_num_prs_components.txt', 'w')
	t.write('tissue\tsample_size\tnum_cafeh_components\tnum_eqtl_components\tavg_num_expressed_genes\n')
	for tissue in tissues:
		samp_size = tissue_to_sample_size[tissue]
		num_comp = tissue_to_num_components[tissue]
		num_eqtls = tissue_to_num_eqtls[tissue]
		avg_num_expressed_genes = tissue_to_number_expressed_genes[tissue]
		t.write(tissue + '\t' + str(samp_size) + '\t' + str(num_comp) + '\t' + str(num_eqtls) + '\t' + str(avg_num_expressed_genes) + '\n')
	t.close()











tissues = []
tissue_to_sample_size = {}
tissue_to_num_components = {}
tissue_to_num_eqtls = {}
tissue_to_number_expressed_genes = {}


f = open(gtex_tissue_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	samp_size = data[1]
	tissues.append(tissue_name)
	tissue_to_sample_size[tissue_name] = samp_size
	tissue_to_num_components[tissue_name] = 0
	tissue_to_num_eqtls[tissue_name] = 0

	composit_tissues = data[3].split(',')
	avg_num_genes = get_avg_number_of_expressed_genes_across_composit_tissues(composit_tissues)
	tissue_to_number_expressed_genes[tissue_name] = avg_num_genes

f.close()
tissues = np.asarray(tissues)


coloc_thresholds = [.1, .3, .5, .7, .9]

for coloc_threshold in coloc_thresholds:
	for tissue in tissues:
		tissue_to_num_components[tissue] = 0
		tissue_to_num_eqtls[tissue] = 0

	for chrom_num in range(1,23):
		component_file = bivariate_cafeh_output_dir + 'mm_v1_coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

		f = open(component_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			num_components = int(data[1])
			num_eqtls = int(data[2])
			tissue_to_num_components[tissue_name] = tissue_to_num_components[tissue_name] + num_components
			tissue_to_num_eqtls[tissue_name] = tissue_to_num_eqtls[tissue_name] + num_eqtls
		f.close()

	# print to output file
	t = open(bivariate_cafeh_output_dir + 'mm_v1_coloc_results_' + trait_name+ '_' + str(coloc_threshold) + '_num_prs_components.txt', 'w')
	t.write('tissue\tsample_size\tnum_cafeh_components\tnum_eqtl_components\tavg_num_expressed_genes\n')
	for tissue in tissues:
		samp_size = tissue_to_sample_size[tissue]
		num_comp = tissue_to_num_components[tissue]
		num_eqtls = tissue_to_num_eqtls[tissue]
		avg_num_expressed_genes = tissue_to_number_expressed_genes[tissue]
		t.write(tissue + '\t' + str(samp_size) + '\t' + str(num_comp) + '\t' + str(num_eqtls) + '\t' + str(avg_num_expressed_genes) + '\n')
	t.close()











tissues = []
tissue_to_sample_size = {}
tissue_to_num_components = {}
tissue_to_num_eqtls = {}
tissue_to_number_expressed_genes = {}


f = open(gtex_tissue_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	samp_size = data[1]
	tissues.append(tissue_name)
	tissue_to_sample_size[tissue_name] = samp_size
	tissue_to_num_components[tissue_name] = 0
	tissue_to_num_eqtls[tissue_name] = 0

	composit_tissues = data[3].split(',')
	avg_num_genes = get_avg_number_of_expressed_genes_across_composit_tissues(composit_tissues)
	tissue_to_number_expressed_genes[tissue_name] = avg_num_genes

f.close()
tissues = np.asarray(tissues)


coloc_thresholds = [.1, .3, .5, .7, .9]

for coloc_threshold in coloc_thresholds:
	for tissue in tissues:
		tissue_to_num_components[tissue] = 0
		tissue_to_num_eqtls[tissue] = 0

	for chrom_num in range(1,23):
		component_file = bivariate_cafeh_output_dir + 'mm_v2_coloc_results_' + trait_name + '_' + str(coloc_threshold) + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

		f = open(component_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			num_components = int(data[1])
			num_eqtls = int(data[2])
			tissue_to_num_components[tissue_name] = tissue_to_num_components[tissue_name] + num_components
			tissue_to_num_eqtls[tissue_name] = tissue_to_num_eqtls[tissue_name] + num_eqtls
		f.close()

	# print to output file
	t = open(bivariate_cafeh_output_dir + 'mm_v2_coloc_results_' + trait_name+ '_' + str(coloc_threshold) + '_num_prs_components.txt', 'w')
	t.write('tissue\tsample_size\tnum_cafeh_components\tnum_eqtl_components\tavg_num_expressed_genes\n')
	for tissue in tissues:
		samp_size = tissue_to_sample_size[tissue]
		num_comp = tissue_to_num_components[tissue]
		num_eqtls = tissue_to_num_eqtls[tissue]
		avg_num_expressed_genes = tissue_to_number_expressed_genes[tissue]
		t.write(tissue + '\t' + str(samp_size) + '\t' + str(num_comp) + '\t' + str(num_eqtls) + '\t' + str(avg_num_expressed_genes) + '\n')
	t.close()