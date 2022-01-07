import numpy as np 
import os
import sys
import pdb


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



trait_name = sys.argv[1]
gtex_tissue_file = sys.argv[2]
processed_bivariate_cafeh_input_dir = sys.argv[3]


gtex_tissues, gtex_tissue_sample_sizes = get_gtex_tissue_names(gtex_tissue_file)


for gtex_tissue in gtex_tissues:
	print(gtex_tissue)
	tissue_merged_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gtex_tissue + '_processed_gene_list.txt'
	t = open(tissue_merged_file, 'w')
	chrom_counter = 0
	for chrom_num in range(1,23):
		input_file = processed_bivariate_cafeh_input_dir + trait_name + '_' + gtex_tissue + '_chr' + str(chrom_num) + '_processed_gene_list.txt'
		f = open(input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			if head_count == 0:
				head_count = head_count + 1
				if chrom_counter == 0:
					t.write(line + '\n')
				continue
			t.write(line + '\n')
		chrom_counter = chrom_counter + 1
		f.close()
	t.close()