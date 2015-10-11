## this is used to plot the corr in cis- multi-linear regression model, between testing dataset (rpkm) and the predicted expression level

import numpy as np
import matplotlib.pyplot as plt



# global variables definition and initialization
num_gene = 0			# TBD

# expression:
gene_list = []			# all genes from the source file
gene_index_map = {}		# re-map those genes into their order (reversed hashing of above)

# information table:
gene_tss = {}			# TSS for all genes (including those pruned genes)
gene_xymt_rep = {}		# map all the X, Y, MT genes

# result table:
corr_rep = {}			# correlation of expected expression level and the real expression level


# information table (colors)
color_table = ['m', '#81b1d2', '#ffed6f', 'r', '#EEEEEE', '#cbcbcb', '#6d904f', 'y', '#E24A33', '#0072B2', '#f0f0f0', '0.40', 'blue', '#fc4f30', '#bfbbd9', '#ccebc4', 'c', '#A60628', '#988ED5', 'g', '#bcbcbc', '#FFB5B8']



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

	###
	gene_tss = {}
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_tss", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]  # don't int() here, as there are XYMT genes
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



	##===================================================== load corr_rep =====================================================
	corr_rep = {}
	file = open("../result_init/para_init_train_cis_corr_dev_ind.txt", 'r')
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
		if gene in gene_xymt_rep:
			continue
		if gene not in corr_rep:
			continue
		corr = corr_rep[gene]
		
		chr = int(gene_tss[gene][0])
		color = color_table[chr-1]
		plt.plot(i, corr, color, marker = 'o', alpha=0.7)

		## add the gene id if the corr is high enough
		if corr >= 0.5:
			print gene,
			print corr


	plt.axis([0, 20000, -1, 1])
	plt.xlabel('Expressed genes (coding and non-coding) from all 22 chromosomes')
	plt.ylabel('Pearson correlation of gene expression level')
	plt.title('Model testing for the multi-linear regression of cis- SNPs (+-1Mb)')
	plt.grid(True)


	plt.show()


