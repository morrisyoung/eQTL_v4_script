## this should be used after the script sample_tissue_preprocess.py
## why above: as the standard to define a NULL gene may be easily updated, while this is not usually the case for eQTL tissues
## this is used to remove NULL genes (according to the standard we defined by ourselves)

## indeed, later on, we can not only remove those NULL genes; we can also remove some non-coding genes, or some specific types of genes


## TODO: this script should be able to quantile normalize all the genes (across all samples)






## the below means "at least \portion of the samples (eQTL) have rpkm value >= \threshold"
threshold = 0.1  ## this rpkm value is used to define a gene is expressed
portion = 0.5  ## 0.9, 0.5


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



	## gene quantile normalization starts from here
