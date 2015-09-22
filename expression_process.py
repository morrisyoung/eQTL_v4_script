## function: further process the expression data, and calculate the principal components for each tissue type
## notes:
##	1. the data is good enough, ranging from "0.11871646" to "346.86022035", which seems good; I will not Gaussian normalize them (maybe not worth)
##	2. whenever I get the new training set, I will need to learn the new parameter initialization

##	3. now we have the expression matrix for each tissue, but we can't use that to initialize the cell state variables (assumed there are 400), as there are very limited sample size (no more than 200) in each of these tissues. so --> I will use all the samples to get the 400 cell state variables (400 principal components here), and use these for all the tissues



import numpy as np
import scipy
import time
from sklearn.decomposition import PCA


num_cellenv = 400



if __name__ == "__main__":

	time_start = time.time()

	##================== get the current training tissues and samples inside (list)
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


	##=================== query the SNP dosage data according to the above list
	




	##=================== calculate the principal components (400) of the training samples
	print "calculating the PCA (400)..."
	expression_matrix = []
	for i in range(len(sample_list)):
		sample = sample_list[i]
		expression_matrix.append(rep_sample[sample])
	X = np.array(expression_matrix)
	pca = PCA(n_components=num_cellenv, copy=True, whiten=False)
	X1 = pca.fit_transform(X)
	#print(pca.explained_variance_ratio_)



	##=============== cellenv to gene
	print "learning the coefficients between PCAs and the expression matrix..."
	U = (np.matrix(X1)).getI() * np.matrix(X)
	#print "the coefficients between them:"
	#print U
	#print "precision comparison:"
	#print X
	#print (np.matrix(X1)) * U
	U = np.squeeze(np.asarray(U))
	U = np.transpose(U)

	## save parameters
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_train_init_cellenv_gene", 'w')
	for i in range(len(U)):
		for j in range(len(U[i])):
			value = U[i][j]
			file.write(str(value) + '\t')
		file.write('\n')
	file.close()



	##=============== snp to cellenv
	print "learning the coefficients between SNP and the PCAs..."



	## save parameters
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_train_init_snp_cellenv", 'w')

	file.close()






	time_end = time.time()
	print "time spent is",
	print time_end - time_start


