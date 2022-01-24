import numpy as np 
import os
import sys
import pdb



def get_snp_positions_from_snp_names(snp_names):
	positions = []
	for snp_name in snp_names:
		snp_info = snp_name.split('_')
		positions.append(int(snp_info[1]))
	return np.asarray(positions)


def load_in_genotype_data(genotype_dosage_file):
	f = open(genotype_dosage_file)
	head_count = 0
	geno_raw = []
	snp_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = np.asarray(data)
			continue
		snp_name = data[0]
		dosages = data[1:]
		snp_names.append(snp_name)
		geno_raw.append(data)
	f.close()
	return np.asarray(header), np.asarray(snp_names), np.asarray(geno_raw)


chrom_num = sys.argv[1]
genotype_dosage_file = sys.argv[2]
cross_tissue_gene_list_file = sys.argv[3]
genotype_reference_panel_dir = sys.argv[4]  # output dir

cis_window = 500000.0

# Load in genotype data
#geno_raw = np.loadtxt(genotype_dosage_file, dtype=str, delimiter='\t')
#header = geno_raw[0,:]
#snp_names = geno_raw[1:,0]
#geno_raw = geno_raw[1:,:]
header, snp_names, geno_raw = load_in_genotype_data(genotype_dosage_file)


# Quick eroror checking
for i,snp_name in enumerate(snp_names):
	if snp_name.endswith('b38') == False:
		print('assumption erororo')
		pdb.set_trace()
	if snp_name != geno_raw[i,0]:
		print('assumptione ororor')
		pdb.set_trace()

# Get snp positions from snp names
snp_positions = get_snp_positions_from_snp_names(snp_names)


f = open(cross_tissue_gene_list_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract fields from this line
	gene_id = data[0]
	line_chrom = data[1]
	# Ignore genes not on desired chromosome
	if line_chrom != 'chr' + chrom_num:
		continue

	tss = int(data[3])

	# Get indices of snps in cis window of gene
	cis_window_snps = np.abs(tss - snp_positions) <= cis_window

	# get genotype data for these cis window snps
	cis_genotype = geno_raw[cis_window_snps,:]
	cis_genotype_to_print = np.vstack((header, cis_genotype))

	# Save to output file
	gene_output_file = genotype_reference_panel_dir + 'genotype_reference_panel_' + gene_id + '.txt'
	np.savetxt(gene_output_file, cis_genotype_to_print, fmt="%s", delimiter='\t')
f.close()