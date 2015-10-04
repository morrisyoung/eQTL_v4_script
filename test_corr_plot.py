## this is used to plot the corr in cis- multi-linear regression model, between testing dataset (rpkm) and the predicted expression level

import numpy as np
import matplotlib.pyplot as plt



# global variables definition and initialization
num_gene = 0			# TBD

# expression:
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
etissue_rep = {}		# mapping all etissues to their indices in filing
gene_chr_map = {}		# we need the chr of genes to diff-color them



# result table:
corr_rep = {}			# correlation of expected expression level and the real expression level


# information table (colors)
color_map = {'1':'m', '2':'#81b1d2', '3':'#ffed6f', '4':'r', '5':'#EEEEEE', '6':'#cbcbcb', '7':'#6d904f', '8':'y', '9':'#E24A33', '10':'#0072B2', '11':'#f0f0f0', '12':'0.40', '13':'blue', '14':'#fc4f30', '15':'#bfbbd9', '16':'#ccebc4', '17':'c', '18':'#A60628', '19':'#988ED5', '20':'g', '21':'#bcbcbc', '22':'#FFB5B8', 'X':'m', 'Y':'r', 'MT':'blue'}



if __name__ == "__main__":


	##===================================================== genes =====================================================
	##=================== gene_list
	file = open("../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized", 'r')
	file.readline()
	file.readline()
	file.readline()
	gene_list = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		gene_list.append(gene)
	file.close()

	###
	gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)
	for i in range(len(gene_list)):
		gene = gene_list[i]
		gene_index_map[gene] = i

	###
	num_gene = len(gene_list)


	##===================================================== plot corr (for each tissue) =====================================================
	## load the etissue_rep first
	etissue_rep = {}
	file = open("../result_predict/etissue_list.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		etissue = line[0]
		etissue_index = int(line[1])
		etissue_rep[etissue] = etissue_index
	file.close()


	## build the chr map
	gene_chr_map = {}
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_tss", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		gene_chr_map[gene] = chr
	file.close()


	for etissue in etissue_rep:
		print etissue
		etissue_index = etissue_rep[etissue]
		file = open("../result_predict/etissue" + str(etissue_index) + ".corr", 'r')
		corr_rep = {}
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			gene = line[0]
			corr = float(line[1])
			corr_rep[gene] = corr
		file.close()


		##===================================================== plot =====================================================
		# with gene_list and corr_rep
		plt.figure(1)
		for i in range(len(gene_list)):
			gene = gene_list[i]
			corr = corr_rep[gene]

			chr = gene_chr_rep[gene]
			color = color_map[chr]
			plt.plot(i, corr, color, marker = 'o', alpha=0.7)

			## add the gene id if the corr is high enough
			if corr >= 0.5:
				print gene,
				print corr


		plt.axis([0, 21000, -1, 1])
		plt.xlabel('Expressed genes (coding and non-coding)')
		plt.ylabel('Pearson correlation of gene expression level')
		plt.title('Testing of modeling precision')
		plt.grid(True)


		plt.show()

