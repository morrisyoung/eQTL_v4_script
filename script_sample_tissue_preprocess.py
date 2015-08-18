## build the qualified sample--tissue rep; to be further used later on
## two criterias: 1. sample has genotype information; and 2. sample is in eQTL (sample size>=60) tissue
## right now I use 60 as the threashold for eQTL tissues

## right after we get the sample list, we should process the rpkm file to only leave these qualified samples (for convenience of up-coming processing and computation)





if __name__ == '__main__':




	'''
	###=================== get the sample attribute (for the tissue property) =======================
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt", 'r')
	file1 = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'w')
	count = 0
	while 1:
		line = (file.readline()).strip()
		count += 1
		if count <= 11:  ## 11
			continue

		if not line:
			break

		line = line.split('\t')
		sample = line[1]
		tissue1 = line[12]
		tissue2 = line[14]

		file1.write(sample + '\t' + tissue1 + '\t' + tissue2 + '\n')

	file.close()
	file1.close()


	## first of all we should remove samples without genotype information
	## #TODO I need to do this, because we have not yet got all the genotype information for all the individuals (overall 550 around).
	# individuals
	individuals = {}
	file = open("GTEx_5M_185_Dec2012.Eigenvectors.txt", 'r')
	file.readline()
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		FID = line[0]
		ID = FID[0:9]
		individuals[ID] = 1
	file.close()


	### check how many tissues and how many samples under each tissue do we have
	file = open("GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file.readline()
	file.readline()
	sample_list = (((file.readline()).strip()).split('\t'))[2:]
	file.close()

	print "there are",
	print len(sample_list),
	print "different samples from the rpkm file."

	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'r')
	rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			#print line
			continue

		sample = line[0]
		tissue = line[2]

		rep[sample] = tissue
	file.close()

	# counting, by the way remove expression data without genotype information
	record = {}
	for sample in sample_list:
		ID = sample[:9]
		if ID not in individuals:
			continue

		tissue = rep[sample]
		if tissue not in record:
			record[tissue] = 1
		else:
			record[tissue] += 1

	# record the tissue-count with size >= 60
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_60", 'w')
	count = 0
	for tissue in record:
		if record[tissue] >= 60:
			count += 1
		#	print tissue,
		#	print record[tissue]
			file.write(tissue + '\t' + str(record[tissue]) + '\n')
	print "# of tissue with sample size >= 60:",
	print count
	file.close()

	# record the tissue-count with size >= 0
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count", 'w')
	count = 0
	for tissue in record:
		count += 1
		#print tissue,
		#print record[tissue]
		file.write(tissue + '\t' + str(record[tissue]) + '\n')
	print "total # of tissue types:",
	print count
	file.close()




	### get all the samples for tissue count >= 60
	## eQTL_tissue
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_60", 'r')
	eQTL_tissue = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		tissue = (line.split('\t'))[0]
		eQTL_tissue[tissue] = []
	file.close()

	## sample_list
	file = open("GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file.readline()
	file.readline()
	sample_list = (((file.readline()).strip()).split('\t'))[2:]
	file.close()

	## sample_tissue_rep
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'r')
	sample_tissue_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			#print line
			continue

		sample = line[0]
		tissue = line[2]

		sample_tissue_rep[sample] = tissue
	file.close()

	## now we have:
	# sample_list
	# sample_tissue_map
	# eQTL_tissue
	for sample in sample_list:
		ID = sample[:9]
		if ID not in individuals:
			continue

		tissue = sample_tissue_rep[sample]
		if tissue in eQTL_tissue:
			eQTL_tissue[tissue].append(sample)
	# sanity check
	#for tissue in rep:
	#	print tissue,
	#	print len(rep[tissue]),
	#	print record[tissue]

	# save the rep
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples", 'w')
	for tissue in eQTL_tissue:
		file.write(tissue + '\t')
		for sample in eQTL_tissue[tissue]:
			file.write(sample + '\t')
		file.write('\n')
	file.close()
	'''




	##============ process the rpkm matrix to get eQTL samples ==============
	## get the sample_rep first
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples", 'r')	
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')[1:]
		for sample in line:
			sample_rep[sample] = 1
	file.close()




	file = open("GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file1 = open("GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed", 'w')
	line = file.readline()
	file1.write(line)
	line = file.readline()
	file1.write(line)


	gene_list_all = []
	index_rep = []



	gene_list = ((file.readline()).strip()).split()

	while 1:
		
