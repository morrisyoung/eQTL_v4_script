## test multi-linear modeling for cis- regulators, with tissue as a variable


import time
import numpy as np


# global variables definition and initialization
num_gene = 0			# TBD
num_tissue = 0			# TBD


# expression:
sample_rep = {}			# hashing all training samples into their rpkm array
sample_list = []		# in-order sample list
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)

sample_tissue_map = {}		# now we need to know the tissue attribute of each sample
tissue_list = []		# we need to know the order of all tissues
tissue_index_map = {}		# reversed hashing of above list


# information table:
gene_xymt_rep = {}		# map all the X, Y, MT genes


# result table:
para_rep = {}
corr_rep = {}			# correlation of expected expression level and the real expression level




if __name__ == "__main__":

	print "working on the multi-linear regression of cis- SNPs for all genes"
	time_start = time.time()



	##===================================================== expression =====================================================
	##================== get the current training tissues and samples inside (list)
	print "get training samples..."
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train", 'r')

	###
	sample_rep = {}
	sample_tissue_map = {}
	tissue_list = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		tissue_list.append(tissue)
		sample_list = line[1:]

		for i in range(len(sample_list)):
			sample = sample_list[i]
			sample_rep[sample] = []
			sample_tissue_map[sample] = tissue

	file.close()

	###
	num_tissue = len(tissue_list)

	###
	tissue_index_map = {}
	for i in range(len(tissue_list)):
		tissue = tissue_list[i]
		tissue_index_map[tissue] = i


	##=================== query the RPKM according to the above list
	print "get expression matrix for these training samples..."
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')

	###
	file.readline()
	file.readline()
	sample_list = ((file.readline()).strip()).split('\t')[2:]
	index_rep = {}
	for i in range(len(sample_list)):
		sample = sample_list[i]
		if sample in sample_rep:
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
				sample = index_rep[i]
				sample_rep[sample].append(rpkm)
	file.close()

	###
	gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
	for i in range(len(gene_list)):
		gene = gene_list[i]
		gene_index_map[gene] = i

	###
	num_gene = len(gene_list)

	###
	sample_list = []
	for sample in sample_rep:
		sample_list.append(sample)

	###
	expression_matrix = []
	for i in range(len(sample_list)):
		sample = sample_list[i]
		expression_matrix.append(sample_rep[sample])
	expression_matrix = np.array(expression_matrix)

	###
	gene_xymt_rep = {}
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_xymt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		gene = line
		gene_xymt_rep[gene] = 1
	file.close()




	##===================================================== regression across all genes =====================================================
	## (we only use the tissue type, with a intercept, to do the regression)
	for i in range(len(gene_list)):

		gene = gene_list[i]
		print gene

		if gene in gene_xymt_rep:  ## we can model the xymt genes simply with tissue types, but to make things comparable, here we don't do that
			continue

		##
		expression_array = []
		for j in range(len(sample_list)):
			rpkm = expression_matrix[j][i]
			expression_array.append(rpkm)
		expression_array = np.array(expression_array)

		##
		genotype_matrix = []
		for j in range(len(sample_list)):
			sample = sample_list[j]
			genotype_matrix.append([])
			'''
			##=== this is for genotype (cis-)
			individual = sample[:9]
			for k in range(start, end+1):
				dosage = snp_dosage_rep[individual][chr-1][k]
				genotype_matrix[j].append(dosage)
			'''
			##=== we need to add the tissue specificity as a variable
			tissue = sample_tissue_map[sample]
			tissue_index = tissue_index_map[tissue]
			list = [0] * num_tissue
			list[tissue_index] = 1
			genotype_matrix[j].extend(list)
			##=== the intercept of the regression
			genotype_matrix[j].append(1)
		genotype_matrix = np.array(genotype_matrix)
		## sample:
		#X = np.array([[1,2,3,1], [2,4,6,1], [3,6,9,1]])  # 1x1 + 2x2 + 3x3 + 1
		#y = np.array([14, 28, 42])
		#y = np.array([15, 29, 44])
		#m = np.linalg.lstsq(X, y)[0]
		m = np.linalg.lstsq(genotype_matrix, expression_array)[0]
		para_rep[gene] = m  ## there is an extra intercept here!!!
		'''
		try:
			m = np.linalg.lstsq(genotype_matrix, expression_array)[0]
			para_rep[gene] = m  ## there is an extra intercept here!!!
		except ValueError:
			print genotype_matrix
			## write the matrix into a file
			file = open("./temp/" + gene + ".txt", 'w')
			for i in range(len(genotype_matrix)):
				for j in range(len(genotype_matrix[i])):
					dosage = genotype_matrix[i][j]
					file.write(str(dosage) + '\t')
				file.write('\n')
			file.close()
		'''



	##===================================================== save all learned parameters =====================================================
	filename = "../result_init/para_init_train_tissuev.txt"
	file = open(filename, 'w')

	for gene in para_rep:
		para_list = para_rep[gene]
		file.write(gene + '\t')
		for i in range(len(para_list)):
			para = para_list[i]
			file.write(str(para) + '\t')
		file.write('\n')
	file.close()







	###
	###
	###
	###
	### part#2
	###
	###
	###
	###
	### note: should use the same tissue_list/tissue_index_map, because the corresponding parameters are for old index






	## get the testing samples first of all
	##===================================================== expression =====================================================
	##================== get the current training tissues and samples inside (list)
	print "get testing samples..."
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test", 'r')
	sample_rep = {}
	sample_tissue_map = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		sample_list = line[1:]

		for i in range(len(sample_list)):
			sample = sample_list[i]
			sample_rep[sample] = []
			sample_tissue_map[sample] = tissue
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
		if sample in sample_rep:
			index_rep[i] = sample

	###
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		rpkm_list = map(lambda x: float(x), line[2:])

		for i in range(len(rpkm_list)):
			rpkm = rpkm_list[i]
			if i in index_rep:
				sample = index_rep[i]
				sample_rep[sample].append(rpkm)
	file.close()


	### this is crucial for building the expression matrix, as we need the order of all samples
	sample_list = []
	for sample in sample_rep:
		sample_list.append(sample)

	###
	expression_matrix = []
	for i in range(len(sample_list)):
		sample = sample_list[i]
		expression_matrix.append(sample_rep[sample])
	expression_matrix = np.array(expression_matrix)



	##===================================================== load all the parameters =====================================================
	filename = "../result_init/para_init_train_tissuev.txt"
	file = open(filename, 'r')
	para_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		para_list = map(lambda x: float(x), line[1:])
		para_rep[gene] = para_list
	file.close()



	##===================================================== testing all genes, get the corr_rep =====================================================
	corr_rep = {}
	for i in range(len(gene_list)):
		gene = gene_list[i]
		print gene

		if gene in gene_xymt_rep:
			continue
		if gene not in para_rep:
			continue

		##
		expression_array_real = []
		for j in range(len(sample_list)):
			rpkm = expression_matrix[j][i]
			expression_array_real.append(rpkm)
		expression_array_real = np.array(expression_array_real)

		##
		expression_array_exp = []
		for j in range(len(sample_list)):
			genotype_array = []
			sample = sample_list[j]
			'''
			##=== we need to add the tissue specificity as a variable
			individual = sample[:9]
			for k in range(start, end+1):
				dosage = snp_dosage_rep[individual][chr-1][k]
				genotype_array.append(dosage)
			'''
			##=== we need to add the tissue specificity as a variable
			tissue = sample_tissue_map[sample]
			tissue_index = tissue_index_map[tissue]
			list = [0] * num_tissue
			list[tissue_index] = 1
			genotype_array.extend(list)
			##=== the intercept of the regression
			genotype_array.append(1)
			##=== expected expression level
			rpkm = 0
			for k in range(len(genotype_array)):
				rpkm += genotype_array[k] * para_rep[gene][k]
			expression_array_exp.append(rpkm)
		expression_array_exp = np.array(expression_array_exp)

		##
		corr = np.corrcoef(expression_array_real, expression_array_exp)[0][1]
		corr_rep[gene] = corr



	##===================================================== save corr_rep =====================================================
	filename = "../result_init/para_init_train_tissuev_corr.txt"
	file = open(filename, 'w')
	for gene in corr_rep:
		file.write(gene + '\t' + str(corr_rep[gene]) + '\n')
	file.close()



	##======== timing
	time_end = time.time()
	print "time spent on this gene is",
	print time_end - time_start


