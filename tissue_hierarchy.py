## to do the tissue hierarchical clustering on all 17 etissues (with both training samples and testing samples)
import numpy as np
import scipy


if __name__ == '__main__':


	'''
	##==================================== == == == == == == == ==
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


	##==================================== == == == == == == == ==
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


	print "number of etissues:",
	print len(tissue_sample_rep)
	print "number of genes:",
	print len(gene_list)
	print "each etissue has the following number of esamples:"
	for tissue in tissue_sample_rep:
		print tissue,
		print len(tissue_sample_rep[tissue])

	file.close()




	##==================================== == == == == == == == ==
	tissue_expression_ave = {}
	for tissue in tissue_sample_rep:
		tissue_expression_ave[tissue] = []
		for i in range(len(gene_list)):
			tissue_expression_ave[tissue].append(0)

		for sample in tissue_sample_rep[tissue]:
			for i in range(len(tissue_sample_rep[tissue][sample])):
				rpkm = tissue_sample_rep[tissue][sample][i]
				tissue_expression_ave[tissue][i] += rpkm

		for i in range(len(gene_list)):
			tissue_expression_ave[tissue][i] = tissue_expression_ave[tissue][i] / len(tissue_sample_rep[tissue])


	print "each etissue has the following number of genes:"
	for tissue in tissue_expression_ave:
		print tissue,
		print len(tissue_expression_ave[tissue])



	##==================================== == == == == == == == ==
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_tissue_mean", 'w')
	for tissue in tissue_expression_ave:
		file.write(tissue + '\t')
		for rpkm in tissue_expression_ave[tissue]:
			file.write(str(rpkm) + '\t')
		file.write('\n')
	print "saving done!"
	'''


	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_tissue_mean", 'r')
	rep = {}
	dimension = 0
	tissue_list = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		expression = map(lambda x: float(x), line[1:])
		rep[tissue] = expression
		dimension = len(rep[tissue])
		tissue_list.append(tissue)
	file.close()

	array = np.ndarray(shape=(len(rep), dimension), dtype=float)
	for i in range(len(tissue_list)):
		tissue = tissue_list[i]
		for j in range(len(rep[tissue])):
			rpkm = rep[tissue][j]
			array[i][j] = rpkm



	Y = scipy.spatial.distance.pdist(array, 'euclidean')
	Z = scipy.cluster.hierarchy.linkage(Y, method='single', metric='euclidean')
	scipy.cluster.hierarchy.dendrogram(Z)
