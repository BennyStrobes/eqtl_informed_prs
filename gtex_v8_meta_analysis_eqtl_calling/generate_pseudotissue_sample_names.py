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
		pseudotissue_name = data[3]
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
	samp_size_dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pseudotissue_name = data[3]
		tissue_name = data[0]
		if pseudotissue_name == 'Remove':
			continue
		if data[2] == 'Remove':
			samp_size_dicti[pseudotissue_name] = data[1]
		if pseudotissue_name not in dicti:
			dicti[pseudotissue_name] = []
		dicti[pseudotissue_name].append(tissue_name)
	f.close()

	pseudotissues = np.sort([*dicti])

	t = open(pseudotissue_info_file,'w')
	t.write('pseudotissue_name\tsample_size\tsample_repeat\tcomposit_tissues\n')

	for pseudotissue in pseudotissues:
		composit_tissues = np.sort(dicti[pseudotissue])
		if pseudotissue not in samp_size_dicti:
			t.write(pseudotissue + '\t320\t')
		else:
			t.write(pseudotissue + '\t' + samp_size_dicti[pseudotissue] + '\t')
		if len(composit_tissues) == 1:
			t.write('False\t')
		else:
			t.write('True\t')
		t.write(','.join(composit_tissues) + '\n')
	t.close()

def get_sample_names_from_composit_tissues(gtex_downsampled_individuals_dir, get_sample_names_from_composit_tissue, composit_tissues):
	sample_names = []
	indi_ids = []
	for tissue in composit_tissues:
		downsample_tissue_file = gtex_downsampled_individuals_dir + 'Downsampled_Individuals_320_' + tissue + '.txt'
		pdb.set_trace()
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

def get_sample_names_from_composit_tissue(gtex_downsampled_individuals_dir, gtex_covariate_dir, tissue, eur_gtex_ind):
	sample_names = []
	indi_ids = []
	downsample_tissue_file = gtex_downsampled_individuals_dir + 'Downsampled_Individuals_320_' + tissue + '.txt'
	if os.path.exists(downsample_tissue_file) == False or tissue == 'Artery_Coronary':
		covariate_file = gtex_covariate_dir + tissue + '.v8.covariates.txt'
		f = open(covariate_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				indi_ids = []
				sample_names = []
				for ele in data:
					if ele in eur_gtex_ind:
						indi_ids.append(ele)
						sample_names.append(tissue + ':' + ele)
				indi_ids = np.asarray(indi_ids)
				sample_names = np.asarray(sample_names)
				continue
			break
		f.close()
	else:
		if tissue == 'Artery_Aorta':
			downsample_tissue_file = gtex_downsampled_individuals_dir + 'Downsampled_Individuals_320_' + tissue + '_v2' + '.txt'
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

def extract_dictionary_list_of_european_ancestry_gtex_individuals(gtex_eur_individuals_file):
	dicti = {}
	f = open(gtex_eur_individuals_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[data[1]] = 1
	f.close()
	return dicti

downsampled_tissue_info_file = sys.argv[1]
gtex_downsampled_individuals_dir = sys.argv[2]
gtex_covariate_dir = sys.argv[3]
pseudotissue_sample_names_dir = sys.argv[4]
gtex_eur_individuals_file = sys.argv[5]


# First generate pseudotissue info file
# One line for each pseudotissue describing its compositional tissues
pseudotissue_info_file = pseudotissue_sample_names_dir + 'pseudotissue_info.txt'
generate_pseudotissue_info_file(downsampled_tissue_info_file, pseudotissue_info_file)


# First generate tissue info file
# One line for each tissue describing it and its assignment to a pseudotissue
tissue_info_file = pseudotissue_sample_names_dir + 'tissue_info.txt'
generate_tissue_info_file(downsampled_tissue_info_file, tissue_info_file)

# Extract dictionary list of european ancestry gtex individuals
eur_gtex_ind = extract_dictionary_list_of_european_ancestry_gtex_individuals(gtex_eur_individuals_file)



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
	sample_names, individual_names = get_sample_names_from_composit_tissue(gtex_downsampled_individuals_dir, gtex_covariate_dir, tissue_name, eur_gtex_ind)
	# Print pseudotissue sample names to output file
	print_arr_to_output(sample_names, pseudotissue_sample_names_dir + tissue_name + '_sample_names.txt')
	print_arr_to_output(individual_names, pseudotissue_sample_names_dir + tissue_name + '_individual_names.txt')
	print_arr_to_output_with_column_of_zeros(individual_names, pseudotissue_sample_names_dir + tissue_name + '_individual_names_plink_ready.txt')

f.close()
