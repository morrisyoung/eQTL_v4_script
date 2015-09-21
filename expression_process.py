## function: further process the expression data, and calculate the principal components for each tissue type

import numpy as np
import scipy



if __name__ == "__main__":

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


		



		break




	file.close()
