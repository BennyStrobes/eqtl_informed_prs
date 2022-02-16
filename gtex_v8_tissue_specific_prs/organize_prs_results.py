import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
from sklearn.linear_model import LinearRegression
import gzip
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import statsmodels.api as sm
from linear_regression_mixtures import LinearRegressionsMixture


def extract_tissue_specific_prs_scores(ukbb_prs_dir, trait_name):
	head_count = 0
	for chrom_num in range(1,23):
		#print(chrom_num)
		file_name = ukbb_prs_dir + trait_name + '_prs_chrom_' + str(chrom_num) + '.txt'
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

def extract_trait_covariates(ukbb_pheno_file1, sample_names):
	sample_dicti = {}
	for sample_name in sample_names:
		sample_dicti[sample_name] = []
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
			trait_names = []
			for i, val in enumerate(data[2:]):
				trait_names.append(val)
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
	trait_cov = np.asarray(arr)
	return trait_cov, trait_names

def extract_covariates_from_cov1(ukbb_pheno_file2, sample_names, cov_name):
	sample_dicti = {}
	for sample_name in sample_names:
		sample_dicti[sample_name] = np.nan
	used_samples = {}
	
	head_count = 0
	f = open(ukbb_pheno_file2)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 95:
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
	reg = LinearRegression(positive=True).fit(prs_train, pheno_train)
	return reg.coef_

def compute_standard_error_from_bootstrap_samples(arr):
	B = len(arr)
	meany = np.mean(arr)
	bootsrap_var = np.sum(np.square(arr-meany))/(B-1)
	return np.sqrt(bootsrap_var)

def learn_prs_weights_wrapper(prs_train, pheno_train, num_boot_strap_samples):
	coef = learn_prs_weights(prs_train, pheno_train)
	
	num_samples = len(pheno_train)
	bootstrap_arr = []
	for bootstrap_sample in np.arange(num_boot_strap_samples):
		bootstrap_sample_indices = np.random.choice(np.arange(num_samples), size=num_samples)
		bootstrap_beta = learn_prs_weights(prs_train[bootstrap_sample_indices,:], pheno_train[bootstrap_sample_indices])
		bootstrap_arr.append(bootstrap_beta)
	bootstrap_arr = np.asarray(bootstrap_arr)
	se = []
	for dimension_num in range(bootstrap_arr.shape[1]):
		dimension_se = compute_standard_error_from_bootstrap_samples(bootstrap_arr[:, dimension_num])
		se.append(dimension_se)
	se = np.asarray(se)
	return coef, se

def compute_relative_r_squared(full_cov_mat, cov_mat, pheno_vec):
	reg = LinearRegression().fit(cov_mat, pheno_vec)
	r_squared = reg.score(cov_mat, pheno_vec)
	reg_full = LinearRegression().fit(full_cov_mat, pheno_vec)
	r_squared_full = reg_full.score(full_cov_mat, pheno_vec)
	return r_squared_full - r_squared

def compute_standard_error_from_jacknife_samples(jacknife_r_squared_arr):
	temp = np.asarray(jacknife_r_squared_arr)
	N_samp = len(temp)
	avg = np.mean(temp)
	vary = ((N_samp-1.0)/(N_samp))*np.sum(np.square(temp-avg))
	return np.sqrt(vary)

def compute_standard_error_from_jacknife_samples_on_subset(jacknife_r_squared_arr, N_tot):
	temp = np.asarray(jacknife_r_squared_arr)
	N_samp = len(temp)
	avg = np.mean(temp)
	vary = ((N_tot-1.0)/(N_tot))*((np.sum(np.square(temp-avg))/N_samp)*N_tot)
	return np.sqrt(vary)	

def compute_relative_r_squared_wrapper(prs_vec_full, pheno_vec_full, cov_mat_full, num_jack_knife_samples):
	X_full = np.asarray(np.hstack((cov_mat_full, np.asmatrix(prs_vec_full).T)))
	X_null = np.asarray(cov_mat_full)
	relative_r_squared = compute_relative_r_squared(X_full, X_null, pheno_vec_full)

	num_samples = len(prs_vec_full)
	jacknife_r_squared_arr = []
	for sample_num in np.random.choice(np.arange(num_samples),size=num_jack_knife_samples, replace=False):
		jacknife_samples = np.delete(np.arange(num_samples), sample_num)
		jacknife_relative_r_squared = compute_relative_r_squared(X_full[jacknife_samples, :], X_null[jacknife_samples,:], pheno_vec_full[jacknife_samples])
		jacknife_r_squared_arr.append(jacknife_relative_r_squared)
	standard_error = compute_standard_error_from_jacknife_samples_on_subset(jacknife_r_squared_arr, num_samples)
	return relative_r_squared, standard_error

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


def extract_technical_covariates(ukbb_pheno_file2,ukbb_pheno_file3, sample_names):
	age = extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'cov_AGE')
	#bmi = extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'body_BMIz')
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
	return cov, cov_names

def get_training_and_testing_indices(num_samples, training_fraction):
	num_training_samples = int(np.round(training_fraction*num_samples))
	training_indices = np.random.choice(num_samples, num_training_samples, replace=False)
	testing_indices = []
	for index in range(num_samples):
		if index not in training_indices:
			testing_indices.append(index)
	testing_indices = np.asarray(testing_indices)
	return training_indices, testing_indices

def print_pca_ve_to_output(pca_ve, output_file):
	t = open(output_file,'w')
	num_pcs = len(pca_ve)
	t.write('pc_num\tvariance_explained\n')
	for pc_num in range(num_pcs):
		t.write('pc' + str(pc_num+1) + '\t' + str(pca_ve[pc_num]) + '\n')
	t.close()

def run_mixture_of_regressions(prs_mat, pheno):
	epsilon = 1e-4
	lam = 0.001
	iterations = 300
	random_restarts = 1
	K=2	
	model = LinearRegressionsMixture(prs_mat, np.expand_dims(pheno, axis=1), K=K)
	model.train(epsilon=epsilon, lam=lam, iterations=iterations, random_restarts=random_restarts, verbose=False)

	return model.probability, model.w, model.pi

def learn_non_linear_prs_weights_wrapper(prs_train, pheno_train, prs_test, pheno_test, cov_test):
	epsilon = 1e-4
	lam = 0.001
	iterations = 2
	random_restarts = 1
	K=2
	model = LinearRegressionsMixture(prs_train, np.expand_dims(pheno_train, axis=1), K=K)
	model.train(epsilon=epsilon, lam=lam, iterations=iterations, random_restarts=random_restarts, verbose=False)

	num_test_samples = prs_test.shape[0]

	predz = []
	posteriors = []

	for test_sample_num in range(num_test_samples):
		y_new, y_posteriors = model.predict(prs_test[test_sample_num,:], posteriors=True)
		predz.append(y_new)
		posteriors.append(y_posteriors)
	predz = np.asarray(predz)
	posteriors = np.asarray(posteriors)

	print(np.sort(posteriors[:,0]))

	#X_full = np.asarray(np.hstack((cov_test, np.asmatrix(predz).T)))
	#X_null = np.asarray(cov_test)
	#relative_r_squared = compute_relative_r_squared(X_full, X_null, pheno_test)

	return predz

def get_tissue_indicees(pseudotissue_file):
	f = open(pseudotissue_file)
	head_count = 0
	counter = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[2] == 'False':
			arr.append(counter)
		counter = counter + 1
	f.close()
	return np.asarray(arr)

ukbb_prs_dir = sys.argv[1]
ukbb_pheno_file1 = sys.argv[2]
ukbb_pheno_file2 = sys.argv[3]
ukbb_pheno_file3 = sys.argv[4]
coloc_threshold = sys.argv[5]
trait_name = sys.argv[6]
analyzed_ukbb_prs_dir = sys.argv[7]
method = sys.argv[8]





output_stem = analyzed_ukbb_prs_dir + trait_name + '_coloc_thresh_' + coloc_threshold + '_' + method + '_'


#pseudotissue_file='/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/pseudotissue_sample_names/pseudotissue_info.txt'
#tissue_indices = get_tissue_indicees(pseudotissue_file)

# Load in data
prs_mat, prs_names, sample_names = extract_tissue_specific_prs_scores(ukbb_prs_dir, trait_name + '_' + coloc_threshold + '_' + method)


#prs_names = prs_names[tissue_indices]
#prs_mat = prs_mat[:, tissue_indices]



# Extract phenotype vector
#pheno =extract_covariates_from_cov2(ukbb_pheno_file2, sample_names, 'blood_WHITE_COUNT_v2')
pheno =extract_covariates_from_cov1(ukbb_pheno_file1, sample_names, trait_name)

# Remove unobserved samples with respect to phenotype of interest
observed = ~np.isnan(pheno)
sample_names = sample_names[observed]
prs_mat = prs_mat[observed,:]
pheno = pheno[observed]
num_samples = len(sample_names)




# Extract covariates
technical_cov, technical_cov_names = extract_technical_covariates(ukbb_pheno_file2,ukbb_pheno_file3, sample_names)
write_matrix_to_output(np.asarray(technical_cov).astype(str), technical_cov_names, sample_names, output_stem + 'technical_covariates.txt')


# Extract blood covariate_names
blood_covariates, blood_covariate_names = extract_blood_trait_covariates(ukbb_pheno_file2, sample_names)
write_matrix_to_output(blood_covariates, blood_covariate_names, sample_names, output_stem + 'blood_covariates.txt')

# Extract trait covariate_names
trait_covariates, trait_covariate_names = extract_trait_covariates(ukbb_pheno_file1, sample_names)
write_matrix_to_output(trait_covariates, trait_covariate_names, sample_names, output_stem + 'trait_covariates.txt')


# Standardize PRS Mat
prs_mat_standardized = StandardScaler().fit_transform(prs_mat)
write_matrix_to_output(prs_mat_standardized.astype(str), np.asarray(prs_names), sample_names, output_stem + 'standardized_prs_scores.txt')
# Remove technical covariates
residual_prs_mat = get_residuals(prs_mat_standardized, technical_cov)
write_matrix_to_output(residual_prs_mat.astype(str), np.asarray(prs_names), sample_names, output_stem + 'standardized_residual_prs_scores.txt')


# Run PCA on prs mat
num_pcs = 10
pca = PCA(n_components=num_pcs)
#principalComponents = pca.fit_transform(residual_prs_mat)
principalComponents = pca.fit_transform(residual_prs_mat)
pca_ve = pca.explained_variance_ratio_
pca_factors = pca.components_


pc_names = []
for pc_num in range(num_pcs):
	pc_names.append('prs_pc' + str(pc_num+1))
write_matrix_to_output(principalComponents.astype(str), np.asarray(pc_names), sample_names, output_stem + 'prs_pca_loadings.txt')
write_matrix_to_output(np.transpose(pca_factors.astype(str)), np.asarray(pc_names), prs_names, output_stem + 'prs_pca_principal_components.txt')
print_pca_ve_to_output(pca_ve, output_stem + 'prs_pca_variance_explained.txt')


# Run mixture of regressions on prs ma
#mor_posteriors, mor_weights, mor_pis  = run_mixture_of_regressions(prs_mat, pheno)
#mor_component_names = np.asarray(['mor_component_1', 'mor_component_2'])
#write_matrix_to_output(mor_posteriors.astype(str), np.asarray(mor_component_names), sample_names, output_stem + 'mixture_of_regressions_posteriors.txt')
#write_matrix_to_output(mor_weights.astype(str), np.asarray(mor_component_names), prs_names, output_stem + 'mixture_of_regressions_weights.txt')
#write_matrix_to_output(np.expand_dims(mor_pis,axis=0).astype(str), np.asarray(mor_component_names), np.asarray(['pi']), output_stem + 'mixture_of_regressions_pi.txt')


# Randomly select training and testing indices
training_indices, testing_indices = get_training_and_testing_indices(num_samples, .2)


# Training data
prs_train = prs_mat[training_indices, :]
prs_train_standardized = prs_mat_standardized[training_indices,:]
pheno_train = pheno[training_indices]
cov_train = technical_cov[training_indices, :]
pc_train = principalComponents[training_indices, :]

# testing data
prs_test = prs_mat[testing_indices, :]
prs_test_standardized = prs_mat_standardized[testing_indices,:]
pheno_test = pheno[testing_indices]
cov_test = technical_cov[testing_indices, :]
pc_test = principalComponents[testing_indices, :]


# Compute correlation between elements of PRS mat
prs_corr_mat = np.corrcoef(np.transpose(prs_mat))
temp = np.hstack((np.asmatrix(prs_names).T, prs_corr_mat.astype(str)))
new_names = []
new_names.append('tissue')
for prs_name in prs_names:
	new_names.append(prs_name)
final = np.vstack((np.asmatrix(new_names), temp))
corr_mat_output_file = output_stem + 'prs_correlation_matrix.txt'
np.savetxt(corr_mat_output_file, final, fmt="%s",delimiter="\t")


# Learn non linear prs weights
#non_linear_joint_prs = learn_non_linear_prs_weights_wrapper(prs_train_standardized, pheno_train, prs_test, pheno_test, cov_test)

# Fit prs weights (to training data)
num_bootstrap_samples = 1000
num_bootstrap_samples = 100
prs_weights, prs_weights_standard_errors = learn_prs_weights_wrapper(prs_train, pheno_train, num_bootstrap_samples)
# Print PRS weights to output
t = open(output_stem + 'prs_weights.txt','w')
t.write('prs_name\tweight\tweight_standard_error\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(prs_weights[index]) + '\t' + str(prs_weights_standard_errors[index]) + '\n')
t.close()

# Fit prs weights (to all samples)
num_bootstrap_samples = 1000
num_bootstrap_samples = 100
prs_weights_all, prs_weights_all_standard_errors = learn_prs_weights_wrapper(prs_mat, pheno, num_bootstrap_samples)
# Print PRS weights to output
t = open(output_stem + 'prs_weights_all_samples.txt','w')
t.write('prs_name\tweight\tweight_standard_error\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(prs_weights_all[index]) + '\t' + str(prs_weights_all_standard_errors[index]) + '\n')
t.close()

# Fit prs weights (to training data)
num_bootstrap_samples = 1000
num_bootstrap_samples = 100
s_prs_weights, s_prs_weights_standard_errors = learn_prs_weights_wrapper(prs_train_standardized, pheno_train, num_bootstrap_samples)
# Print PRS weights to output
t = open(output_stem + 'standardized_prs_weights.txt','w')
t.write('prs_name\tweight\tweight_standard_error\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(s_prs_weights[index]) + '\t' + str(s_prs_weights_standard_errors[index]) + '\n')
t.close()

# Fit prs weights (to all samples)
num_bootstrap_samples = 1000
num_bootstrap_samples = 100
s_prs_weights_all, s_prs_weights_all_standard_errors = learn_prs_weights_wrapper(prs_mat_standardized, pheno, num_bootstrap_samples)
# Print PRS weights to output
t = open(output_stem + 'standardized_prs_weights_all_samples.txt','w')
t.write('prs_name\tweight\tweight_standard_error\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(s_prs_weights_all[index]) + '\t' + str(s_prs_weights_all_standard_errors[index]) + '\n')
t.close()



print('done')

# Relative R-squared
num_jack_knife_samples = 200
num_jack_knife_samples = 2
relative_r_squared_arr = []
relative_r_squared_std_err_arr = []
for prs_index, prs_name in enumerate(prs_names):
	relative_r_squared, relative_r_squared_std_err = compute_relative_r_squared_wrapper(prs_test[:, prs_index], pheno_test, cov_test, num_jack_knife_samples)
	relative_r_squared_arr.append(relative_r_squared)
	relative_r_squared_std_err_arr.append(relative_r_squared_std_err)
# joint model
joint_prs = np.dot(prs_test, prs_weights)
relative_r_squared, relative_r_squared_std_err = compute_relative_r_squared_wrapper(joint_prs, pheno_test, cov_test, num_jack_knife_samples)
relative_r_squared_arr.append(relative_r_squared)
relative_r_squared_std_err_arr.append(relative_r_squared_std_err)

# Print relative r-squared to output
t = open(output_stem + 'relative_r_squared.txt','w')
t.write('prs_name\trelative_r_squared\trelative_r_squared_standard_error\n')
for index,prs_name in enumerate(prs_names):
	t.write(prs_name + '\t' + str(relative_r_squared_arr[index]) + '\t' + str(relative_r_squared_std_err_arr[index]) + '\n')
t.write('joint_prs\t' + str(relative_r_squared_arr[-1]) + '\t' + str(relative_r_squared_std_err_arr[-1]) + '\n')

#non_linear_r_squared, non_linear_r_squared_std_err = compute_relative_r_squared_wrapper(non_linear_joint_prs, pheno_test, cov_test, num_jack_knife_samples)
#t.write('non_linear_joint_prs\t' + str(non_linear_r_squared) + '\t' + str(non_linear_r_squared_std_err) + '\n')
t.close()


# Relative r-squared with PCs
t = open(output_stem + 'relative_r_squared_with_prs_pcs.txt','w')
t.write('prs_name\trelative_r_squared\trelative_r_squared_standard_error\n')
for pc_num in range(pc_test.shape[1]):
	print('PC' + str(pc_num))
	print(np.corrcoef(pc_test[:, pc_num], pheno_test)[0,1])
	relative_r_squared, relative_r_squared_std_err = compute_relative_r_squared_wrapper(pc_test[:, pc_num], pheno_test, cov_test, num_jack_knife_samples)
	t.write('prs_pc' + str(pc_num + 1) + '\t' + str(relative_r_squared) + '\t' + str(relative_r_squared_std_err) + '\n')
t.close()



