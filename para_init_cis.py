## this is used to initialzie the cis- parameters (actually learn these parameters from multi-linear regression)

## the learned parameters have no tissue specificity (that's why this is only the initialization)

## the results are saved in "../result_init/para_init_train_cis.txt", each line is for one gene (several SNPs are in the cis- region of this gene)

## we can later on use the same script to test the prediction precision, and draw the plot for all the 20603 genes





# global variables definition and initialization
num_snp = 0;			# TBD
num_cellenv = 400;		# Specified
num_gene = 0;			# TBD
num_etissue = 0;		# TBD
num_batch = 0;			# TBD
num_batch_hidden = 100;		# Specified
num_individual = 0;		# TBD


## TODO: not yet filled in the above numbers


# genotype:
snp_name_list = []
snp_pos_list = []


# expression:
eQTL_tissue_rep = {}		# hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
eQTL_samples = {}		# hashing all eQTL samples to their tissues
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
etissue_list = []		# eTissues in order
etissue_index_map = {}		# re-map those etissues into their order (reversed hashing of above)
esample_tissue_rep = {}		# esample lists of all etissues

# information table:
gene_tss = {}			# TSS for all genes (including those pruned genes)
gene_xymt_rep = {}		# map all the X, Y, MT genes


# batch:
batch_individual = {}
batch_sample = {}








def genotype_read(snp_dosage_list, individual):
	## fill in the snp_dosage_list (read from file)
	for i in range(22):
		file = open("../genotype_185_dosage_matrix_qc/chr" + str(i+1) + "/SNP_dosage_" + individual + ".txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			dosage = float(line)
			snp_dosage_list[i].append(dosage)
		file.close()
	return




if __name__ == "__main__":


	print "working on the multi-linear regression of cis- SNPs for all genes"





	## read in all the relevant data, and prepare for prediction/showcase of the learned parameters
	##===================================================== genotype =====================================================
	print "snp_name_list and snp_pos_list..."
	snp_name_list = []
	snp_pos_list = []
	for i in range(22):
		file = open("../genotype_185_dosage_matrix_qc/chr" + str(i+1) + "/SNP_info.txt", 'r')
		snp_name_list.append([])
		snp_pos_list.append([])
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')
			snp = line[0]
			pos = int(line[1])
			snp_name_list[-1].append(snp)
			snp_pos_list[-1].append(pos)
		file.close()

	'''
	snp_dosage_list = []
	for i in range(22):
		snp_dosage_list.append([])
	genotype_read(snp_dosage_list, "GTEX-QDVN")
	'''

	##===================================================== expression =====================================================
	print "esample_tissue_rep..."
	esample_tissue_rep = {}
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		etissue = line[0]
		sample_list = line[1:]
		esample_tissue_rep[etissue] = sample_list
	file.close()

	print "eQTL_samples..."
	eQTL_samples = {}
	for etissue in esample_tissue_rep:
		for i in range(len(esample_tissue_rep[etissue])):
			esample = esample_tissue_rep[etissue][i]
			eQTL_samples[esample] = etissue

	print "etissue_list and etissue_index_map..."
	etissue_list = []
	etissue_index_map = {}
	file = open("../result/etissue_list.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break
		
		line = line.split('\t')
		etissue = line[0]
		etissue_index = int(line[1]) - 1
		etissue_list.append(etissue)
		etissue_index_map[etissue] = etissue_index
	file.close()

	print "eQTL_tissue_rep, gene_list and gene_index_map..."
	eQTL_tissue_rep = {}
	gene_list = []
	gene_index_map = {}
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')
	file.readline()
	file.readline()
	line = ((file.readline()).strip()).split('\t')[2:]
	rep = {}
	for i in range(len(line)):
		sample = line[i]
		if sample in eQTL_samples:
			etissue = eQTL_samples[sample]
			if etissue not in eQTL_tissue_rep:
				eQTL_tissue_rep[etissue] = {}
			eQTL_tissue_rep[etissue][sample] = []
			rep[i] = sample
	index = 0
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		expr_list = line[2:]
		gene_list.append(gene)
		gene_index_map[gene] = index
		index += 1

		for i in range(len(expr_list)):
			if i in rep:
				esample = rep[i]
				etissue = eQTL_samples[esample]
				rpkm = float(expr_list[i])
				eQTL_tissue_rep[etissue][esample].append(rpkm)
	file.close()
	
	print "gene_tss..."
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

	print "gene_xymt_rep..."
	gene_xymt_rep = {}
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_xymt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		gene = line
		gene_xymt_rep[gene] = 1
	file.close()


	##===================================================== batch =====================================================
	print "batch_individual..."
	batch_individual = {}
	file = open("../batch_var_individual.txt", 'r')
	file.readline()
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		individual = line[0]
		batch_list = map(lambda x: float(x), line[1:])
		batch_individual[individual] = batch_list
	file.close()

	print "batch_sample..."
	batch_sample = {}
	file = open("../batch_var_sample.txt", 'r')
	file.readline()
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		sample = line[0]
		batch_list = map(lambda x: float(x), line[1:])
		batch_sample[sample] = batch_list
	file.close()


	##===================================================== parameters =====================================================
	'''
	## parameter containers:
	vector<vector<float *>> para_cis_gene;
	vector<float *> para_snp_cellenv;
	vector<vector<float *>> para_cellenv_gene;
	vector<float *> para_batch_batch_hidden;
	vector<float *> para_batch_hidden_gene;

	# information table:
	unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position in the snp vector)
	'''

	print "para_cis_gene..."
	para_cis_gene = []
	for i in range(len(etissue_list)):
		para_cis_gene.append([])
	for i in range(len(para_cis_gene)):
		etissue = etissue_list[i]
		print i+1,
		print etissue
		index = i+1
		file = open("../result/para_cis_gene/etissue" + str(index) + ".txt", 'r')
		# each line is a gene
		while 1:
			line = file.readline()
			if not line:
				break

			line = line.strip()
			if line == '':
				para_cis_gene[i].append([])
			else:
				line = line.split('\t')
				line = map(lambda x: float(x), line)
				para_cis_gene[i].append(line)
		file.close()

	print "para_snp_cellenv..."
	file = open("../result/para_snp_cellenv.txt", 'r')


	'''
	para_snp_cellenv = []
	'''
	

	



