import numpy as np 
import os
import sys
import pdb



def extract_ordered_genes_from_gene_list(cafeh_gene_list):
	f = open(cafeh_gene_list)
	genes = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes.append(data[0])
	return np.asarray(genes)

def association_file_to_bed_file(cafeh_hg38_file, temp_hg38_bed):
	f = open(cafeh_hg38_file)
	t = open(temp_hg38_bed,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 5:
			print('assumption eroror')
			pdb.set_trace()
		snp_id = data[2]
		snp_info = snp_id.split('_')
		chrom_num = snp_info[0]
		pos = snp_info[1]
		pos2 = int(pos) + 1
		t.write(chrom_num + '\t' + pos + '\t' + str(pos2) + '\n')
	f.close()
	t.close()

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg38ToHg19.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)


def recreate_cafeh_association_file_in_hg19_format(hg38_association_file, hg19_association_file, liftover_output_file, liftover_missing_file, hg38_to_hg19_mapping_file):
	missing = {}
	f = open(liftover_missing_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		missing[data[0] + '_' + data[1]] = 1
	f.close()
	f = open(hg38_association_file)
	t = open(hg19_association_file,'w')
	t_mapping = open(hg38_to_hg19_mapping_file,'w')
	g = open(liftover_output_file)
	t_mapping.write('hg38_snp_id\thg19_snp_id\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		snp_id = data[2]
		snp_info = snp_id.split('_')
		snp_key = snp_info[0] + '_' + snp_info[1]
		if snp_key in missing:
			t_mapping.write(snp_id + '\tNA\n')
			continue
		hg19_info = next(g).rstrip().split('\t')
		new_snp_id = hg19_info[0] + '_' + hg19_info[1] + '_' + snp_info[2] + '_' + snp_info[3]
		t.write(','.join(data[:2]) + ',' + new_snp_id + ',' + ','.join(data[3:]) + '\n')
		t_mapping.write(snp_id + '\t' + new_snp_id + '\n')
	f.close()
	g.close()
	t.close()
	t_mapping.close()

def get_num_lines(file_name):
	counter = 0
	f = open(file_name)
	for line in f:
		counter = counter + 1
	f.close()
	return counter

def error_checking(liftover_output_file, hg19_association_file):
	num_lines_lift = get_num_lines(liftover_output_file)
	num_lines_association = get_num_lines(hg19_association_file) - 1
	if num_lines_association != num_lines_lift:
		print('assumption erororo')


def get_genes_on_chromosome(cafeh_gene_list, chrom_num):
    f = open(cafeh_gene_list)
    counter = 0
    head_count = 0
    genes = []
    gene_dicti = {}
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        counter = counter + 1
        if head_count == 0:
            head_count = head_count + 1
            continue
        gene_name = data[0]
        line_chrom_num = data[1]
        if line_chrom_num != 'chr' + chrom_num:
            continue
        genes.append(gene_name)
        gene_dicti[gene_name] = []
    f.close()
    return np.asarray(genes)



cafeh_gene_list = sys.argv[1]
processed_gtex_associations_dir = sys.argv[2]
liftover_directory = sys.argv[3]
chrom_num = sys.argv[4]


# Extract ordered list of genes on this chromosome
ordered_genes = get_genes_on_chromosome(cafeh_gene_list, chrom_num)


# Loop through genes
for gene_index, gene in enumerate(ordered_genes):
	# Cafeh file in hg38 format
	cafeh_hg38_file = processed_gtex_associations_dir + gene + '_associations.csv'
	# First convert hg38  association to temporary bed format
	temp_hg38_bed = processed_gtex_associations_dir + gene + '_temp_bed_hg38.txt'
	association_file_to_bed_file(cafeh_hg38_file, temp_hg38_bed)

	# Run liftover
	liftover_output_file = processed_gtex_associations_dir + gene + '_liftover_output.txt'
	liftover_missing_file = processed_gtex_associations_dir + gene + '_liftover_missing.txt'
	run_liftover(temp_hg38_bed, liftover_output_file, liftover_missing_file, liftover_directory)


	# Cafeh association file in hg19 format
	cafeh_hg19_file = processed_gtex_associations_dir + gene + '_hg19_associations.csv'
	hg38_to_hg19_mapping_file = processed_gtex_associations_dir + gene + '_hg38_to_hg19_snp_mapping.txt'
	recreate_cafeh_association_file_in_hg19_format(cafeh_hg38_file, cafeh_hg19_file, liftover_output_file, liftover_missing_file, hg38_to_hg19_mapping_file)

	# Quick error checking
	error_checking(liftover_output_file, cafeh_hg19_file)


	# Remove temporary files
	os.system('rm ' + liftover_output_file)
	os.system('rm ' + liftover_missing_file)
	os.system('rm ' + temp_hg38_bed)

