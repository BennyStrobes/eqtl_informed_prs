import numpy as np 
import os
import sys
import pdb
from sklearn.linear_model import LinearRegression
import gzip
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def extract_tissue_specific_prs_scores(ukbb_prs_dir, trait_name, threshold):
	head_count = 0
	for chrom_num in range(1,23):
		#print(chrom_num)
		file_name = ukbb_prs_dir + trait_name + '_cafeh_prs_beta_threshold_' + threshold + '_chrom_' + str(chrom_num) + '.txt'
		raw_data = np.loadtxt(file_name, dtype=str,delimiter='\t')
		tissue_names = raw_data[0,1:]
		sample_names = raw_data[1:,0]
		prs = raw_data[1:,1:].astype(float)
		#print(np.std(prs))
		#print(np.mean(prs))
		if head_count == 0:
			head_count = head_count + 1
			prs_mat = np.copy(prs)
			tissue_names_arr = np.copy(tissue_names)
			sample_names_arr = np.copy(sample_names)
		else:
			prs_mat = prs_mat + prs
			if np.array_equal(tissue_names, tissue_names_arr) == False or np.array_equal(sample_names, sample_names_arr) == False:
				print('assumption eroorr')
				pdb.set_trace()
	return prs_mat, tissue_names_arr, sample_names_arr

def extract_blood_trait_covariates(ukbb_pheno_file2, sample_names):
	sample_dicti = {}
	for sample_name in sample_names:
		sample_dicti[sample_name] = []
	used_samples = {}
	head_count = 0
	f = open(ukbb_pheno_file2)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 45:
			print('assumption eroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			trait_names = []
			for i, val in enumerate(data[2:]):
				trait_names.append(val.split('_v2')[0])
			continue
		sample_name = data[0]
		if sample_name not in sample_dicti:
			continue
		sample_dicti[sample_name] = np.asarray(data[2:])
		used_samples[sample_name] = 1
	f.close()
	arr = []
	for sample_name in sample_names:
		arr.append(sample_dicti[sample_name])
	blood_cov = np.asarray(arr)
	return blood_cov, trait_names

def extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, cov_name):
	sample_dicti = {}
	for sample_name in sample_names:
		sample_dicti[sample_name] = np.nan
	used_samples = {}
	
	head_count = 0
	f = open(ukbb_pheno_file2)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 45:
			print('assumption eroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			pheno_index = np.where(np.asarray(data) == cov_name)[0][0]
			continue	
		sample_name = data[0]
		if sample_name not in sample_dicti:
			continue
		used_samples[sample_name] = 1
		if data[pheno_index] == 'NA':
			continue
		sample_dicti[sample_name] = float(data[pheno_index])
	f.close()

	if len(used_samples) != len(sample_dicti):
		print('assumption eroorr')
		pdb.set_trace()

	arr = []
	for sample_name in sample_names:
		arr.append(sample_dicti[sample_name])

	if len(arr) != len(sample_dicti):
		print('assumption erorror')
		pdb.set_trace()

	return np.asarray(arr)



def extract_pheno(ukbb_pheno_file1, sample_names, pheno_name):
	sample_dicti = {}
	for sample_name in sample_names:
		sample_dicti[sample_name] = np.nan
	used_samples = {}

	head_count = 0
	f = open(ukbb_pheno_file1)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 95:
			print('assumption eroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			pheno_index = np.where(np.asarray(data) == pheno_name)[0][0]
			continue
		sample_name = data[0]
		if sample_name not in sample_dicti:
			continue
		used_samples[sample_name] = 1
		if data[pheno_index] == 'NA':
			continue
		sample_dicti[sample_name] = float(data[pheno_index])
	f.close()

	if len(used_samples) != len(sample_dicti):
		print('assumption eroorr')
		pdb.set_trace()

	arr = []
	for sample_name in sample_names:
		arr.append(sample_dicti[sample_name])

	if len(arr) != len(sample_dicti):
		print('assumption erorror')
		pdb.set_trace()

	return np.asarray(arr)

def extract_genotype_pcs_from_cov3(ukbb_pheno_file3, sample_names, num_pcs):
	sample_dicti = {}
	for sample_name in sample_names:
		sample_dicti[sample_name] = np.nan
	used_samples = {}

	f = gzip.open(ukbb_pheno_file3)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if len(data) != 112:
			print('assumtpion eroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_name = data[0]
		if sample_name not in sample_dicti:
			continue
		pcs = np.asarray(data[72:(72+num_pcs)]).astype(float)
		sample_dicti[sample_name] = pcs
		used_samples[sample_name] = 1
	f.close()

	if len(used_samples) != len(sample_dicti):
		print('assumption eroorr')
		pdb.set_trace()

	arr = []
	for sample_name in sample_names:
		arr.append(sample_dicti[sample_name])

	return np.asarray(arr)


def learn_prs_weights(prs_train, pheno_train):
	reg = LinearRegression(positive=True).fit(prs_mat, pheno)
	return reg.coef_

def compute_relative_r_squared(prs_vec, pheno_vec, cov_mat):
	full_cov_mat = np.hstack((cov_mat, np.asmatrix(prs_vec).T))

	reg = LinearRegression().fit(np.asarray(cov_mat), pheno_vec)
	r_squared = reg.score(np.asarray(cov_mat), pheno_vec)
	
	reg_full = LinearRegression().fit(np.asarray(full_cov_mat), pheno_vec)
	r_squared_full = reg_full.score(np.asarray(full_cov_mat), pheno_vec)
	return r_squared_full - r_squared

def run_pca_on_prs_mat(prs_mat):
	pdb.set_trace()

def temper(pc_vec, blood_covariates, blood_covariate_names):
	for index, blood_covariate_name in enumerate(blood_covariate_names):
		print(blood_covariate_name)
		blood_covariate_vec = blood_covariates[:, index]
		observed = blood_covariate_vec != 'NA'
		blood_covariate_vec2 = (blood_covariate_vec[observed]).astype(float)

		print(np.corrcoef(blood_covariate_vec2, pc_vec[observed])[0,1])

def write_matrix_to_output(data_mat, column_names, row_names, output_file):
	t = open(output_file,'w')
	t.write('sample_name\t' + '\t'.join(column_names) + '\n')
	for row_index, row_name in enumerate(row_names):
		t.write(row_name + '\t' + '\t'.join(data_mat[row_index,:]) + '\n')
	t.close()


def get_residuals(prs_mat, cov):
	prs_resid_mat = np.copy(prs_mat)
	for prs_num in range(prs_mat.shape[1]):
		reg = LinearRegression().fit(np.asarray(cov), prs_mat[:, prs_num])
		prs_resid_mat[:, prs_num] = prs_mat[:, prs_num] - reg.predict(np.asarray(cov))
	prs_resid_mat = StandardScaler().fit_transform(prs_resid_mat)

	return prs_resid_mat


def make_variable_categorical(var):
	num_samples = len(var)
	num_categories = len(np.unique(var))
	unique_categories = np.unique(var)
	mapping = {}
	for i, category in enumerate(unique_categories):
		mapping[category] = i
	new_data = np.zeros((num_samples, num_categories))
	for sample_index, sample_var in enumerate(var):
		new_data[sample_index,mapping[sample_var]] =1 
	return new_data[:,:(num_categories-1)]


ukbb_prs_dir = sys.argv[1]
ukbb_pheno_file1 = sys.argv[2]
ukbb_pheno_file2 = sys.argv[3]
ukbb_pheno_file3 = sys.argv[4]
thresh = sys.argv[5]


trait_name = 'blood_white_count'


output_stem = ukbb_prs_dir + trait_name + '_'




# Load in data
prs_mat, prs_names, sample_names = extract_tissue_specific_prs_scores(ukbb_prs_dir, trait_name, thresh)




# Extract phenotype vector
pheno =extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'blood_WHITE_COUNT_v2')




# Remove unobserved samples
observed = ~np.isnan(pheno)
sample_names = sample_names[observed]
prs_mat = prs_mat[observed,:]

# Extract covariates
age = extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'cov_AGE')
bmi = extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'body_BMIz')
sex = extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'cov_SEX')
access_center = extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'cov_ASSESS_CENTER')
access_center_categorical = make_variable_categorical(access_center)

genotype_pcs = extract_genotype_pcs_from_cov3(ukbb_pheno_file3, sample_names, 10)
cov = np.hstack((np.asmatrix(age).T, np.asmatrix(sex).T, genotype_pcs, access_center_categorical))
cov_names = []
cov_names.append('age')
cov_names.append('sex')
for pc_num in range(10):
	cov_names.append('genotype_pc' + str(pc_num+1))
for assess_center_num in range(access_center_categorical.shape[1]):
	cov_names.append('assess_center' + str(assess_center_num+1))
cov_names = np.asarray(cov_names)

# Extract blood covariate_names
blood_covariates, blood_covariate_names = extract_blood_trait_covariates(ukbb_pheno_file2, sample_names)
write_matrix_to_output(blood_covariates, blood_covariate_names, sample_names, output_stem + 'blood_covariates.txt')

# Standardize PRS Mat
prs_mat = StandardScaler().fit_transform(prs_mat)

residual_prs_mat = get_residuals(prs_mat, cov)

# Run PCA on prs mat
num_pcs = 5
pca = PCA(n_components=num_pcs)
principalComponents = pca.fit_transform(residual_prs_mat)
pc_names = []
for pc_num in range(num_pcs):
	pc_names.append('prs_pc' + str(pc_num+1))
write_matrix_to_output(principalComponents.astype(str), np.asarray(pc_names), sample_names, output_stem + 'prs_pcs.txt')




pheno = pheno[observed]
num_samples = len(sample_names)

# Randomly select training and testing indices
num_training_samples = int(np.round(.05*num_samples))
training_indices = np.random.choice(num_samples, num_training_samples, replace=False)
testing_indices = []
for index in range(num_samples):
	if index not in training_indices:
		testing_indices.append(index)
testing_indices = np.asarray(testing_indices)

# Training data
prs_train = prs_mat[training_indices, :]
pheno_train = pheno[training_indices]
cov_train = cov[training_indices, :]
pc_train = principalComponents[training_indices, :]

# testing data
prs_test = prs_mat[testing_indices, :]
pheno_test = pheno[testing_indices]
cov_test = cov[testing_indices, :]
pc_test = principalComponents[testing_indices, :]



for index in range(cov.shape[1]):
	print(cov_names[index])
	print(np.corrcoef(np.asarray(cov)[:,index], principalComponents[:,1])[0,1])
for index in range(cov.shape[1]):
	print(cov_names[index])
	print(np.corrcoef(np.asarray(cov)[:,index], principalComponents[:,2])[0,1])
for index in range(cov.shape[1]):
	print(cov_names[index])
	print(np.corrcoef(np.asarray(cov)[:,index], principalComponents[:,0])[0,1])

prs_corr_mat = np.corrcoef(np.transpose(prs_mat))
temp = np.hstack((np.asmatrix(prs_names).T, prs_corr_mat.astype(str)))
new_names = []
new_names.append('tissue')
for prs_name in prs_names:
	new_names.append(prs_name)
final = np.vstack((np.asmatrix(new_names), temp))
corr_mat_output_file = output_stem + 'prs_correlation_matrix.txt'
np.savetxt(corr_mat_output_file, final, fmt="%s",delimiter="\t")

print(prs_corr_mat)
# Fit prs weights
prs_weights = learn_prs_weights(prs_train, pheno_train)
print(prs_weights)
pdb.set_trace()

# Print relative r-squared to output
t = open(output_stem + 'prs_weights.txt','w')
t.write('prs_name\tweight\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(prs_weights[index]) + '\n')
t.close()







relative_r_squared_arr = []
for prs_index, prs_name in enumerate(prs_names):
	print(prs_name)
	print(np.corrcoef(prs_test[:,prs_index], pheno_test)[0,1])
	relative_r_squared = compute_relative_r_squared(prs_test[:, prs_index], pheno_test, cov_test)
	print(relative_r_squared)
	relative_r_squared_arr.append(relative_r_squared)

print('joint prs')
joint_prs = np.dot(prs_test, prs_weights)
print(np.corrcoef(joint_prs, pheno_test)[0,1])
relative_r_squared = compute_relative_r_squared(joint_prs, pheno_test, cov_test)
print(relative_r_squared)
relative_r_squared_arr.append(relative_r_squared)

# PCS
for pc_num in range(pc_test.shape[1]):
	print('PC' + str(pc_num))
	print(np.corrcoef(pc_test[:, pc_num], pheno_test)[0,1])
	relative_r_squared = compute_relative_r_squared(pc_test[:, pc_num], pheno_test, cov_test)
	print(relative_r_squared)


# Print relative r-squared to output
t = open(output_stem + 'relative_r_squared.txt','w')
t.write('prs_name\trelative_r_squared\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(relative_r_squared_arr[index]) + '\n')
t.write('joint_prs\t' + str(relative_r_squared_arr[-1]) + '\n')
t.close()




