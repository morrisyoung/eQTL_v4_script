## to do the tissue hierarchical clustering on all 17 etissues (with both training samples and testing samples)


if __name__ == '__main__':


	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type", 'r')
	sample_tissue_map = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		if len(line) < 3:
			continue
		sample = line[0]
		tissue = line[2]
		sample_tissue_map[sample] = tissue
	file.close()



	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')
	file.readline()
	file.readline()
	line = (file.readline()).strip()
	line = (line.split('\t'))[2:]
	sample_list = line
	tissue_list = []
	tissue_sample_rep = {}
	for sample in sample_list:
		tissue = sample_tissue_map[sample]
		tissue_list.append(tissue)
		if tissue in tissue_sample_rep:
			tissue_sample_rep[tissue][sample] = []
		else:
			tissue_sample_rep[tissue] = {}
			tissue_sample_rep[tissue][sample] = []
	gene_list = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		line = line[2:]
		line = map(lambda x: float(x), line)

		gene_list.append(gene)
		for i in range(len(line)):
			rpkm = line[i]
			sample = sample_list[i]
			tissue = tissue_list[i]
			tissue_sample_rep[tissue][sample].append(rpkm)


	print len(tissue_sample_rep)
	print len(gene_list)
	for tissue in tissue_sample_rep:
		print tissue,
		print len(tissue_sample_rep[tissue])

	file.close()

