import numpy as np 
import os
import sys
import pdb






gtex_tissue_file = sys.argv[1]
trait_name = sys.argv[2]
bivariate_cafeh_output_dir = sys.argv[3]



tissues = []
tissue_to_sample_size = {}
tissue_to_num_components = {}
tissue_to_num_eqtls = {}


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
f.close()
tissues = np.asarray(tissues)


for chrom_num in range(1,23):
	component_file = bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_num_prs_components_chrom_' + str(chrom_num) + '.txt'

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
t = open(bivariate_cafeh_output_dir + 'cafeh_results_' + trait_name + '_num_prs_components.txt', 'w')
t.write('tissue\tsample_size\tnum_cafeh_components\tnum_cafeh_eqtl_components\n')
for tissue in tissues:
	samp_size = tissue_to_sample_size[tissue]
	num_comp = tissue_to_num_components[tissue]
	num_eqtls = tissue_to_num_eqtls[tissue]
	t.write(tissue + '\t' + str(samp_size) + '\t' + str(num_comp) + '\t' + str(num_eqtls) + '\n')
t.close()

