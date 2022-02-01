import numpy as np 
import os
import sys
import pdb



def make_expression_data_matrix_eqtl_ready_in_one_chromosome(expression_file, matrix_eqtl_expression_file, matrix_eqtl_gene_location_file, chrom_num):
	chrom_string = 'chr' + str(chrom_num)

	f = open(expression_file)
	t = open(matrix_eqtl_expression_file,'w')
	t2 = open(matrix_eqtl_gene_location_file,'w')

	t2.write('geneid\tchr\ts1\ts2\n')

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('id\t' + '\t'.join(data[4:]) + '\n')
			continue
		line_chrom_num = data[1]
		# ignore genes not on current chromosome
		if line_chrom_num != chrom_string:
			continue
		gene_id = data[0]
		gene_tss = data[2]
		t.write(gene_id + '\t' + '\t'.join(data[4:]) + '\n')
		t2.write(gene_id + '\t' + line_chrom_num + '\t' + gene_tss + '\t' + gene_tss + '\n')
	t.close()
	t2.close()
	f.close()



expression_file = sys.argv[1]
matrix_eqtl_expression_output_root = sys.argv[2]
matrix_eqtl_gene_location_output_root = sys.argv[3]


for chrom_num in range(1,23):
	matrix_eqtl_expression_file = matrix_eqtl_expression_output_root + 'chr' + str(chrom_num) + '.txt'
	matrix_eqtl_gene_location_file = matrix_eqtl_gene_location_output_root + 'chr' + str(chrom_num) + '.txt'

	make_expression_data_matrix_eqtl_ready_in_one_chromosome(expression_file, matrix_eqtl_expression_file, matrix_eqtl_gene_location_file, chrom_num)