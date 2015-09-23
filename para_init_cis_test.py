## this is used to test the init parameters of the cis- xxx
## I will calculate the predicted rpkm for all samples, and compare that with the true rpkm, and draw the comparison plot
## the learned parameters have no tissue specificity (that's why this is only the initialization)


import time
import numpy as np
import matplotlib.pyplot as plt


# global variables definition and initialization
num_gene = 0			# TBD
num_individual = 0		# TBD

# genotype:
snp_pos_list = []		# the position of all SNPs
snp_dosage_rep = {}		# load all the dosage (for cluster)

# expression:
sample_rep = {}			# hashing all training samples into their rpkm array
sample_list = []		# in-order sample list
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)

# information table:
gene_tss = {}			# TSS for all genes (including those pruned genes)
gene_xymt_rep = {}		# map all the X, Y, MT genes
gene_cis_index = {}		# mapping the gene to cis snp indices (start position and end position in the snp vector)

# result table:
para_rep = {}			# the learned cis- parameters for all the non-xymt genes
corr_rep = {}			# correlation of expected expression level and the real expression level





def snp_dosage_load():
	global snp_dosage_rep

	for individual in snp_dosage_rep:
		for i in range(22):
			snp_dosage_rep[individual].append([])
			file = open("../genotype_185_dosage_matrix_qc/chr" + str(i+1) + "/SNP_dosage_" + individual + ".txt", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break
 
				dosage = float(line)
				snp_dosage_rep[individual][i].append(dosage)
			file.close()
	return




if __name__ == "__main__":

	print "testing the multi-linear parameters..."
	time_start = time.time()


	##===================================================== genotype =====================================================
	print "snp_dosage_rep..."
	file = open("../list_individual.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line
		snp_dosage_rep[individual] = []
	file.close()

	num_individual = len(snp_dosage_rep)
	snp_dosage_load()


	## snp_pos_list
	print "snp_pos_list..."
	snp_pos_list = []
	for i in range(22):
		snp_pos_list.append([])
		file = open("../genotype_185_dosage_matrix_qc/chr" + str(i+1) + "/SNP_info.txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')
			snp = line[0]
			pos = int(line[1])
			snp_pos_list[i].append(pos)
		file.close()



	##===================================================== expression =====================================================
	##================== get the current training tissues and samples inside (list)
	print "get testing samples..."
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test", 'r')
	sample_rep = {}
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
	gene_tss = {}
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_tss", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		tss = int(line[2])
		gene_tss[gene] = (chr, tss)
	file.close()

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




	##===================================================== cis- region definition =====================================================
	gene_cis_index = {}
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_cis_range", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		start = int(line[1])
		end = int(line[2])
		gene_cis_index[gene] = (start, end)
	file.close()



	##===================================================== load all the parameters =====================================================
	file = open("../result_init/para_init_train_cis.txt", 'r')
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



	##===================================================== testing all genes =====================================================
	corr_rep = {}
	for i in range(len(gene_list)):

		gene = gene_list[i]
		print gene

		if gene in gene_xymt_rep:
			continue

		if gene not in para_rep:
			continue

		chr = int(gene_tss[gene][0])
		start = gene_cis_index[gene][0]
		end = gene_cis_index[gene][1]

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
			individual = sample[:9]
			for k in range(start, end+1):
				dosage = snp_dosage_rep[individual][chr-1][k]
				genotype_array.append(dosage)
			genotype_array.append(1)  # we need the intercept
			## expected expression level
			rpkm = 0
			for k in range(len(genotype_array)):
				rpkm += genotype_array[k] * para_rep[gene][k]
			expression_array_exp.append(rpkm)
		expression_array_exp = np.array(expression_array_exp)

		##
		corr = np.corrcoef(expression_array_real, expression_array_exp)[0][1]
		corr_rep[gene] = corr




	##===================================================== plot =====================================================
	# with gene_list and corr_rep
	plt.figure(1)
	for i in range(len(gene_list)):
		gene = gene_list[i]
		if gene in gene_xymt_rep:
			continue
		if gene not in corr_rep:
			continue
		corr = corr_rep[gene]
		plt.plot(i, corr, 'bo')

	plt.show()




	##======== timing
	time_end = time.time()
	print "time spent on this gene is",
	print time_end - time_start



