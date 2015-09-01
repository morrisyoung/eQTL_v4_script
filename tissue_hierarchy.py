## to do the tissue hierarchical clustering on all 17 etissues (with both training samples and testing samples)
import numpy as np
import scipy
from scipy import spatial
from scipy import cluster
from matplotlib import pyplot as plt
import dendropy
from dendropy.calculate import treemeasure
import math


pair_list_rep = {}
tissue_list = []




def func_dfs(Tree):
	global pair_list_rep

	ID = Tree.get_id()
	dist = Tree.dist

	if Tree.is_leaf():
		return [ID]

	else:
		Tree_left = Tree.get_left()
		Tree_right=  Tree.get_right()
		list_left = func_dfs(Tree_left)
		list_left = list_left[:]
		list_right = func_dfs(Tree_right)
		list_right = list_right[:]

		pair_list_rep[ID] = (dist, list_left, list_right)

		list_return = []
		for element in list_left:
			list_return.append(element)
		for element in list_right:
			list_return.append(element)

		return list_return






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
	#file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r') # TODO
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene", 'r')
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
	#file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_tissue_mean", 'w') # TODO
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_tissue_mean", 'w')
	for tissue in tissue_expression_ave:
		file.write(tissue + '\t')
		for rpkm in tissue_expression_ave[tissue]:
			file.write(str(rpkm) + '\t')
		file.write('\n')
	print "saving done!"


	'''







	##================================== from mean expression of genes to hierarchical clustering ==================================
	#file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized_tissue_mean", 'r') # TODO
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_tissue_mean", 'r')
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
		for j in range(dimension):
			rpkm = rep[tissue][j]
			array[i][j] = rpkm


	##============ hierarchical clustering ============
	Y = spatial.distance.pdist(array, 'euclidean')
	Z = cluster.hierarchy.linkage(Y, method='single', metric='euclidean')	
	##=========== get the pairwise distance matrix ===========
	Tree = cluster.hierarchy.to_tree(Z)
	func_dfs(Tree)

	#file = open("../tissue_hierarchy_normalized.txt", 'w') # TODO
	file = open("../tissue_hierarchy_unnormalized.txt", 'w')
	for node in pair_list_rep:
		dist = pair_list_rep[node][0]
		list_left = pair_list_rep[node][1]
		list_right = pair_list_rep[node][2]

		for element1 in list_left:
			for element2 in list_right:
				print tissue_list[element1],
				print " <--> ",
				print tissue_list[element2],
				print " : ",
				print dist
				file.write(tissue_list[element1] + '\t' + tissue_list[element2] + '\t' + str(dist) + '\n')
	file.close()
	##========================================================


	R = cluster.hierarchy.dendrogram(Z, orientation='left', labels=tissue_list, color_threshold=0)


	plt.xlabel('distance')
	plt.show()

