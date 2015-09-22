## function: further process the expression data, and calculate the principal components for each tissue type
## notes:
##	1. the data is good enough, ranging from "0.11871646" to "346.86022035", which seems good; I will not Gaussian normalize them (maybe not worth)
##	2. whenever I get the new training set, I will need to learn the new parameter initialization

##	3. now we have the expression matrix for each tissue, but we can't use that to initialize the cell state variables (assumed there are 400), as there are very limited sample size (no more than 200) in each of these tissues. so --> I will use all the samples to get the 400 cell state variables (400 principal components here), and use these for all the tissues

##	4. currently I use the following inversion method: M_{observe} = M_{pca} * M_{coefficient} --> M_{pca}^{-1} * M_{observe} = M_{coefficient}



import numpy as np
import scipy
import time
from sklearn.decomposition import PCA


num_cellenv = 400



if __name__ == "__main__":

	time_start = time.time()

	##================== get the current training tissues and samples inside (list)
	print "get training samples..."
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train", 'r')
	rep_sample = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		sample_list = line[1:]


		for i in range(len(sample_list)):
			sample = sample_list[i]
			rep_sample[sample] = []

	file.close()


	##=================== query the RPKM according to the above list
	print "get expression matrix for these training samples..."
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')

	file.readline()
	file.readline()
	sample_list = ((file.readline()).strip()).split('\t')[2:]
	rep_index = {}
	for i in range(len(sample_list)):
		sample = sample_list[i]
		if sample in rep_sample:
			rep_index[i] = sample

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = (line.split('\t'))[2:]
		rpkm_list = map(lambda x: float(x), line)

		for i in range(len(rpkm_list)):
			rpkm = rpkm_list[i]
			if i in rep_index:
				sample = rep_index[i]
				rep_sample[sample].append(rpkm)
	file.close()
	sample_list = []
	for sample in rep_sample:
		sample_list.append(sample)
	expression_matrix = []
	for i in range(len(sample_list)):
		sample = sample_list[i]
		expression_matrix.append(rep_sample[sample])
	expression_matrix = np.array(expression_matrix)



	''' don't need this any more
	##=================== query the SNP dosage data according to the above list
	print "get the SNP dosage for these training samples..."
	rep_snp = {}
	for i in range(len(sample_list)):
		sample = sample_list[i]
		individual = sample[:9]
		if individual in rep_snp:
			continue
		else:
			rep_snp[individual] = []

	for individual in rep_snp:
		## looping the 22 chr
		print individual
		for i in range(22):
			chr = i+1
			file = open("../genotype_185_dosage_matrix_qc/chr" + str(chr) + "/SNP_dosage_" + individual + ".txt", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				dosage = float(line)
				rep_snp[individual].append(dosage)
			file.close()
	snp_matrix = []
	for i in range(len(sample_list)):
		sample = sample_list[i]
		individual = sample[:9]
		snp_list = rep_snp[individual][:]
		snp_matrix.append(snp_list)
	snp_matrix = np.array(snp_matrix)
	'''
	



	##=================== calculate the principal components (400) of the training samples
	print "calculating the PCA (400, or num_cellenv)..."
	pca = PCA(n_components=num_cellenv, copy=True, whiten=False)
	Y = pca.fit_transform(expression_matrix)
	#print(pca.explained_variance_ratio_)






	## the folloiwng are for the coefficients between cellenv and gene, and snp and cellenv
	## this will take "22131.8078818" seconds to finish (including the above code, which takes only a little time)
	'''
	##=============== cellenv to gene
	print "learning the coefficients between PCAs and the expression matrix..."
	U = (np.matrix(Y)).getI() * np.matrix(expression_matrix)
	U = np.squeeze(np.asarray(U))
	U = np.transpose(U)

	## save parameters
#	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_train_init_cellenv_gene", 'w')
	for i in range(len(U)):
		for j in range(len(U[i])):
			value = U[i][j]
			file.write(str(value) + '\t')
		file.write('\n')
	file.close()


	##=============== snp to cellenv
	print "learning the coefficients between SNP and the PCAs..."
	U = (np.matrix(snp_matrix)).getI() * (np.matrix(Y))
	U = np.squeeze(np.asarray(U))
	U = np.transpose(U)

	# test: shuould be num_cellenv, num_snp
	print len(U)
	print len(U[0])

	## save parameters
#	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_train_init_snp_cellenv", 'w')
	for i in range(len(U)):
		for j in range(len(U[i])):
			value = U[i][j]
			file.write(str(value) + '\t')
		file.write('\n')
	file.close()
	'''






	time_end = time.time()
	print "time spent is",
	print time_end - time_start


