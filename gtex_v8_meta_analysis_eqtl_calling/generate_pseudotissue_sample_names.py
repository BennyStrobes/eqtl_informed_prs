import numpy as np 
import os
import sys
import pdb


def generate_tissue_info_file(downsampled_tissue_info_file, tissue_info_file):
	f = open(downsampled_tissue_info_file)
	head_count = 0
	dicti = {}
	t = open(tissue_info_file,'w')
	t.write('tissue_name\tsample_size\tpseudotissue_name\n')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pseudotissue_name = data[2]
		tissue_name = data[0]
		if pseudotissue_name == 'Remove':
			continue
		samp_size = int(data[1])
		t.write(tissue_name + '\t' + str(samp_size) + '\t' + pseudotissue_name + '\n')
	f.close()
	t.close()


def generate_pseudotissue_info_file(downsampled_tissue_info_file, pseudotissue_info_file):
	f = open(downsampled_tissue_info_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pseudotissue_name = data[2]
		tissue_name = data[0]
		if pseudotissue_name == 'Remove':
			continue
		if pseudotissue_name not in dicti:
			dicti[pseudotissue_name] = []
		dicti[pseudotissue_name].append(tissue_name)
	f.close()

	pseudotissues = np.sort([*dicti])

	t = open(pseudotissue_info_file,'w')
	t.write('pseudotissue_name\tsample_size\tsample_repeat\tcomposit_tissues\n')

	for pseudotissue in pseudotissues:
		composit_tissues = np.sort(dicti[pseudotissue])
		t.write(pseudotissue + '\t320\t')
		if len(composit_tissues) == 1:
			t.write('False\t')
		else:
			t.write('True\t')
		t.write(','.join(composit_tissues) + '\n')
	t.close()

def get_sample_names_from_composit_tissues(gtex_downsampled_individuals_dir, composit_tissues):
	sample_names = []
	indi_ids = []
	for tissue in composit_tissues:
		downsample_tissue_file = gtex_downsampled_individuals_dir + 'Downsampled_Individuals_320_' + tissue + '.txt'
		f = open(downsample_tissue_file)
		for line in f:
			indi_id = line.rstrip()
			sample_name = tissue + ':' + indi_id
			sample_names.append(sample_name)
			indi_ids.append(indi_id)
		f.close()
	sample_names = np.asarray(sample_names)
	indi_ids = np.asarray(indi_ids)
	if len(sample_names) != 320:
		print('assumption eroror')
		pdb.set_trace()
	if len(indi_ids) != 320:
		print('assumption erororo')
		pdb.set_trace()
	return sample_names, indi_ids

def get_sample_names_from_composit_tissue(gtex_downsampled_individuals_dir, tissue):
	sample_names = []
	indi_ids = []
	downsample_tissue_file = gtex_downsampled_individuals_dir + 'Downsampled_Individuals_320_' + tissue + '.txt'
	f = open(downsample_tissue_file)
	for line in f:
		indi_id = line.rstrip()
		sample_name = tissue + ':' + indi_id
		sample_names.append(sample_name)
		indi_ids.append(indi_id)
	f.close()
	sample_names = np.asarray(sample_names)
	indi_ids = np.asarray(indi_ids)
	return sample_names, indi_ids


def print_arr_to_output(arr, output_file):
	t = open(output_file,'w')
	for ele in arr:
		t.write(ele + '\n')
	t.close()

def print_arr_to_output_with_column_of_zeros(arr, output_file):
	t = open(output_file,'w')
	for ele in arr:
		t.write('0\t' + ele + '\n')
	t.close()

def print_sample_repeat_structure(individual_names, output_file):
	dicti = {}
	for counter, individual_name in enumerate(np.sort(np.unique(individual_names))):
		dicti[individual_name] = counter
	t = open(output_file,'w')
	for individual_name in individual_names:
		t.write(str(dicti[individual_name]) + '\n')
	t.close()

downsampled_tissue_info_file = sys.argv[1]
gtex_downsampled_individuals_dir = sys.argv[2]
pseudotissue_sample_names_dir = sys.argv[3]



# First generate pseudotissue info file
# One line for each pseudotissue describing its compositional tissues
pseudotissue_info_file = pseudotissue_sample_names_dir + 'pseudotissue_info.txt'
generate_pseudotissue_info_file(downsampled_tissue_info_file, pseudotissue_info_file)

# First generate tissue info file
# One line for each tissue describing it and its assignment to a pseudotissue
tissue_info_file = pseudotissue_sample_names_dir + 'tissue_info.txt'
generate_tissue_info_file(downsampled_tissue_info_file, tissue_info_file)



# Loop through tissues and print ordered sample names for that tissue
head_count = 0
f = open(tissue_info_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# get relevent fields from line
	tissue_name = data[0]

	# Extract ordered list of sample names from this pseudotissue
	sample_names, individual_names = get_sample_names_from_composit_tissue(gtex_downsampled_individuals_dir, tissue_name)
	# Print pseudotissue sample names to output file
	print_arr_to_output(sample_names, pseudotissue_sample_names_dir + tissue_name + '_sample_names.txt')
	print_arr_to_output(individual_names, pseudotissue_sample_names_dir + tissue_name + '_individual_names.txt')
	print_arr_to_output_with_column_of_zeros(individual_names, pseudotissue_sample_names_dir + tissue_name + '_individual_names_plink_ready.txt')

f.close()
