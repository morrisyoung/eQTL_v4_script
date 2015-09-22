## this is used to initialzie the cis- parameters (actually learn these parameters from multi-linear regression)

## the learned parameters have no tissue specificity (that's why this is only the initialization)

## the results are saved in "../result_init/para_init_train_cis.txt", each line is for one gene (several SNPs are in the cis- region of this gene)

## we can later on use the same script to test the prediction precision, and draw the plot for all the 20603 genes



import time
import numpy as np





# global variables definition and initialization
num_gene = 0;			# TBD
num_individual = 0;		# TBD

# genotype:
snp_pos_list = []		# the position of all SNPs
snp_dosage_rep = {}		# hash all individuals to their dosage list (a list of 22 lists)

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
para_rep = {}




def snp_dosage_load():
	global snp_dosage_rep

	for inidividual in snp_dosage_rep:
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

	print "working on the multi-linear regression of cis- SNPs for all genes"
	time_start = time.time()

	##===================================================== genotype =====================================================
	## snp_dosage_rep
	print "snp_dosage_rep..."
	snp_dosage_rep = {}
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


	###
	### this takes 149.607182026 seconds
	###



	##===================================================== expression =====================================================
	##================== get the current training tissues and samples inside (list)
	print "get training samples..."
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train", 'r')
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
	# gene_cis_index
	for i in range(len(gene_list)):
		gene = gene_list[i]

		if gene in gene_xymt_rep:
			continue

		chr = int(gene_tss[gene][0])
		tss = gene_tss[gene][1]

		flag1 = 0
		flag2 = 0
		start = 0
		end = 0
		for j in range(len(snp_pos_list[chr-1])):
			if flag1 == 0:
				if (snp_pos_list[chr-1][j] - tss >= -1000000) and (snp_pos_list[chr-1][j] - tss <= 1000000):
					start = j
					flag1 = 1
			if flag1 == 1 and flag2 == 0:
				if (snp_pos_list[chr-1][j] - tss >= -1000000) and (snp_pos_list[chr-1][j] - tss <= 1000000):
					end = j
				else:
					flag2 = 1;
			if (flag1 == 1) and (flag2 == 1):
				break

		gene_cis_index[gene] = (start, end)

	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_cis_range", 'r')
	for gene in gene_cis_index:
		start = gene_cis_index[gene][0]
		end = gene_cis_index[gene][1]
		file.write(gene + '\t' + str(start) + '\t' + str(end) + '\n')
	file.close()



	'''
	##===================================================== regression across all genes =====================================================
	for i in range(len(gene_list)):
		gene = gene_list[i]

		if gene in gene_xymt_rep:
			continue

		para_rep[gene] = []  # store the learned parameters here
	'''

		
		








	time_end = time.time()
	print "time spent is",
	print time_end - time_start




