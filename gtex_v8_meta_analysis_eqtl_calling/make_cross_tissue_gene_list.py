import numpy as np 
import os
import sys
import pdb


def extract_tissue_names(pseudotissue_info_file):
	arr = []
	f = open(pseudotissue_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[0]
		arr.append(tissue_name)
	return np.asarray(arr)


pseudotissue_info_file = sys.argv[1]
pseudotissue_expression_dir = sys.argv[2]


# Extract tissue names
tissue_names = extract_tissue_names(pseudotissue_info_file)

# Create dictionary with word for every gene used along with mapping to vector of tissues
# also gene to tss mapping
gene_to_tissues = {}
gene_to_tss = {}


for tissue in tissue_names:
	print(tissue)
	for chrom_num in range(1,23):
		gene_loc_file = pseudotissue_expression_dir + tissue + '_gene_location_matrix_eqtl_ready_chr' + str(chrom_num) + '.txt'
		f = open(gene_loc_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			gene_id = data[0]
			chrom_num = data[1]
			pos = data[2]
			if gene_id not in gene_to_tissues:
				gene_to_tissues[gene_id] = []
			gene_to_tissues[gene_id].append(tissue)
			if gene_id not in gene_to_tss:
				gene_to_tss[gene_id] = (chrom_num, pos)
			else:
				if gene_to_tss[gene_id][0] != chrom_num or gene_to_tss[gene_id][1] != pos:
					print('assumption eorororr')
					pdb.set_trace()
		f.close()

genes = np.sort([*gene_to_tissues])
t = open(pseudotissue_expression_dir + 'cross_tissue_gene_list.txt','w')
t.write('gene_name\tchrom_num\tactive_tissues\thg38_tss\n')
for gene_id in genes:
	chrom_num = gene_to_tss[gene_id][0]
	tss = gene_to_tss[gene_id][1]
	tissue_arr = np.sort(np.asarray(gene_to_tissues[gene_id]))
	t.write(gene_id + '\t' + chrom_num + '\t' + ';'.join(tissue_arr) + '\t' + tss + '\n')
t.close()
