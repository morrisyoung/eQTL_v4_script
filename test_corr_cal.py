## calculate the Pearson correlation on testing dataset


import time
import numpy as np
import matplotlib.pyplot as plt



# global variables definition and initialization
num_gene = 0			# TBD
num_etissue = 0			# TBD


# expression:
esample_rep = {}		# hashing all testing samples into their etissue
eQTL_rep = {}			# etissue: (esample1:[xx,xx.xx], esample2:[xx,xx,xx])
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
esample_tissue_rep = {}		# ordered samples in each tissue
etissue_rep = {}		# read from results


# result table:
corr_rep = {}			# correlation of expected expression level and the real expression level




if __name__ == "__main__":

	print "testing the multi-linear parameters..."
	time_start = time.time()


	##===================================================== expression (real) =====================================================
	##================== get the current training tissues and samples inside (list)
	print "get testing samples..."
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test", 'r')
	sample_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		etissue = line[0]
		esample_list = line[1:]
		
		eQTL_rep[etissue] = {}

		for i in range(len(esample_list)):
			esample = esample_list[i]
			esample_rep[esample] = etissue
			eQTL_rep[etissue][esample] = []

	file.close()


	##=================== query the RPKM according to the above list
	print "get expression matrix for these testing samples..."
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')

	###
	file.readline()
	file.readline()
	sample_list = ((file.readline()).strip()).split('\t')[2:]
	index_rep = {}
	for i in range(len(sample_list)):
		sample = sample_list[i]
		if sample in esample_rep:
			index_rep[i] = sample

	###
	gene_list = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		gene_list.append(gene)
		rpkm_list = map(lambda x: float(x), line[2:])

		for i in range(len(rpkm_list)):
			rpkm = rpkm_list[i]
			if i in index_rep:
				esample = index_rep[i]
				etissue = esample_rep[esample]
				eQTL_rep[etissue][esample].append(rpkm)
	file.close()

	###
	gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
	for i in range(len(gene_list)):
		gene = gene_list[i]
		gene_index_map[gene] = i

	###
	num_gene = len(gene_list)

	### build the sample list for each tissue; this is crucial for building the expression matrix, as we need the order of all samples
	esample_tissue_rep = {}
	for etissue in eQTL_rep:
		esample_tissue_rep[etissue] = []
		for esample in eQTL_rep[etissue]:
			esample_tissue_rep[etissue].append(esample)



	##===================================================== corr (for each tissue) =====================================================
	##===== we will work on each tissue
	## get the etissue_rep first
	file = open("../result_predict/etissue_list.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		etissue = line[0]
		etissue_index = int(line[1])
		etissue_rep[etissue] = etissue_index
	file.close()


	num_etissue = len(etissue_rep)


	for etissue in etissue_rep:
		print etissue
		etissue_index = etissue_rep[etissue]
		## read and save the predicted gene expression
		rep = {}
		file = open("../result_predict/etissue" + str(etissue_index) + ".txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			sample = line[0]
			rpkm_list = map(lambda x: float(x), line[1:])
			rep[sample] = rpkm_list
		file.close()

		## calculate the Pearson corr
		corr_rep = {}
		for i in range(len(gene_list)):
			gene = gene_list[i]

			##
			expression_array_real = []
			for j in range(len(esample_tissue_rep[etissue])):
				esample = esample_tissue_rep[etissue][j]
				rpkm = eQTL_rep[etissue][esample][i]
				expression_array_real.append(rpkm)
			expression_array_real = np.array(expression_array_real)

			##
			expression_array_exp = []
			for j in range(len(esample_tissue_rep[etissue])):
				esample = esample_tissue_rep[etissue][j]
				rpkm = rep[esample][i]
				expression_array_exp.append(rpkm)
			expression_array_exp = np.array(expression_array_exp)

			##
			corr = np.corrcoef(expression_array_real, expression_array_exp)[0][1]
			corr_rep[gene] = corr


		file = open("../result_predict/etissue" + str(etissue_index) + ".corr", 'w')
		for gene in corr_rep:
			file.write(gene + '\t' + str(corr_rep[gene]) + '\n')
		file.close()



	##======== timing
	time_end = time.time()
	print "time spent on this gene is",
	print time_end - time_start

