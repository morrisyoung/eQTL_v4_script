## function: further process the expression data, and calculate the principal components for each tissue type
## notes:
##	1. the data is good enough, ranging from "0.11871646" to "346.86022035", which seems good; I will not Gaussian normalize them (maybe not worth)
##	2. whenever I get the new training set, I will need to learn the new parameter initialization


import numpy as np
import scipy



if __name__ == "__main__":


	##================== get the current training tissues and samples inside (list)
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train", 'r')
	rep_sample = {}
	rep_tissue_sample = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		sample_list = line[1:]

		rep_tissue_sample[tissue] = {}

		for i in range(len(sample_list)):
			sample = sample_list[i]
			rep_sample[sample] = tissue
			rep_tissue_sample[tissue][sample] = []

	file.close()



	##=================== query the RPKM according to the above list
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')

	file.readline()
	file.readline()
	sample_list = ((file.readline()).strip()).split('\t')[2:]


	rpkm_matrix = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = (line.split('\t'))[2:]
		line = map(lambda x: float(x), line)

		## we need to pick out the training data (simply drop the testing data)
		rpkm_matrix.append(line)

		## we are all working on the designated training dataset currently
	

		break

	file.close()




