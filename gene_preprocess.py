## this should be used after the script sample_tissue_preprocess.py, as this will work on the gene dimension
## why above: as the standard to define a NULL gene may be easily updated, while this is not usually the case for eQTL tissues
## this script is used to:
## 1. remove NULL genes (according to the standard we defined by ourselves; indeed, we can also remove some non-coding genes, or some specific types of genes)
## 2. quantile normalize all the genes across all selected samples

import numpy as np




## NULL gene removal standard:
## the below means "at least \portion of the samples (eQTL) have rpkm value >= \threshold"
threshold = 0.1  ## this rpkm value is used to define a gene is expressed
portion = 0.5  ## 0.9, 0.5


## quantile-normalization:
## the standard method



if __name__ == '__main__':



	'''
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_1_sample", 'r')
	file1 = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene", 'w')

	line = file.readline()
	file1.write(line)
	line = file.readline()
	file1.write(line)
	line = file.readline()
	file1.write(line)

	count = 0  # total number of genes
	count1 = 0  # selected number of genes
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		description = line[1]
		expression = map(lambda x: float(x), line[2:])

		total = len(expression)
		num = 0
		for rpkm in expression:
			if rpkm >= threshold:
				num += 1
		ratio = num * 1.0 / total

		if ratio >= portion:
			count1 += 1
			file1.write(gene + '\t' + description + '\t')
			for rpkm in expression:
				file1.write(str(rpkm) + '\t')
			file1.write('\n')

		count += 1

	file.close()
	file1.close()

	print "there are totally",
	print count,
	print "genes, and",
	print count1,
	print "of them will be included into our eQTL analysis (have at least",
	print portion,
	print "of all samples that have",
	print threshold,
	print "rpkm value)."
	'''


	##========================== gene quantile normalization starts from here ==========================
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene", 'r')
	file1 = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'w')

	line = file.readline()
	file1.write(line)
	line = file.readline()
	file1.write(line)
	line = file.readline()
	file1.write(line)


	gene_list = []
	description_list = []
	rpkm_rep = {}  # hashing gene to L1
	L2 = []

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		## sorting; indexing (additive filling L2)
		line = line.split('\t')
		gene = line[0]
		description = line[1]
		gene_list.append(gene)
		description_list.append(description)

		expression = map(lambda x: float(x), line[2:])
		expression = np.array(expression)

		sort = np.argsort(expression)
		## get the ordering list
		sort1 = []
		for i in range(len(sort)):
			sort1.append(0)
		for i in range(len(sort)):
			index = sort[i]
			sort1[index] = i

		rpkm_rep[gene] = sort1

		if len(L2) == 0:
			for pos in sort:
				rpkm = expression[pos]
				L2.append(rpkm)
		else:
			for i in range(len(sort)):
				pos = sort[i]
				rpkm = expression[pos]
				L2[i] += rpkm
	file.close()

	length = len(gene_list)
	for i in range(len(L2)):
		L2[i] = L2[i] * 1.0 / length


	for i in range(len(gene_list)):
		## two lists:
		## L1 (value as the re-mapped positions of all original elements; each gene has one such list)
		## L2 (containing the normalized/averaged value for each index/position)
		gene = gene_list[i]
		description = description_list[i]
		L1 = rpkm_rep[gene]

		file1.write(gene + '\t' + description + '\t')

		for index in L1:
			value = L2[index]
			file1.write(str(value) + '\t')

		file1.write('\n')

	file1.close()
