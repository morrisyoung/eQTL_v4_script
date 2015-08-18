## some scripts for testing something
## this is really the try.py script




if __name__ == '__main__':

	'''
	##==== get the gene positions in chromosome#22 ====
	file = open("gencode.v18.genes.patched_contigs.gtf", 'r')
	file1 = open("gencode.v18.genes.patched_contigs.gtf_chr22", 'w')
	count = 0

	while 1:
		line = (file.readline()).strip()
		count += 1
	
		if count <= 5:
			continue

		if not line:
			break

		line = line.split('\t')

		if line[0] == '22' and line[2] == 'transcript':
			for element in line:
				file1.write(element + '\t')
			file1.write('\n')

	file.close()
	file1.close()
	'''



	'''
	##==== get the gene (1171) rpkm for chromsome#22 ====
	file = open("GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file1 = open("gencode.v18.genes.patched_contigs.gtf_chr22", 'r')
	file2 = open("GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_chr22", 'w')


	## get all the genes
	gene_rep = {}
	while 1:
		line = (file1.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[8][9: 26]
		gene_rep[gene] = 1

	file1.close()
	

	count = 0

	number = 0
	while 1:
		line = (file.readline()).strip()
		count += 1

		if count <= 2:
			continue

		if not line:
			break

		if count == 3:
			file2.write(line + '\n')
			continue

		gene = line[0:17]
		if gene in gene_rep:
			file2.write(line + '\n')
			number += 1

	file.close()
	file2.close()
	'''









	'''

	## testing whether all the recorded samples (eQTL samples) have their genotype information
	# individuals
	individuals = {}
	file = open("./GTEx_5M_185_Dec2012.Eigenvectors.txt", 'r')
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


	count = 0
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		for i in range(1, len(line)):
			sample = line[i]
			ID = sample[:9]
			if ID not in individuals:
				print sample
				count += 1
	print count

	'''



	''' test how many samples we have from tissue_sample-count list
	file = open("phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_60", 'r')
	total = 0
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		count = int(line[1])
		total = total + count

	file.close()
	print total
	'''
