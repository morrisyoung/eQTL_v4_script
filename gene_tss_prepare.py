## process the annotation file to get the tss for all genes first


if __name__ == '__main__':

	## get the chr and tss for all genes (including X, Y and MT genes)
	file = open("../gencode.v18.genes.patched_contigs.gtf", 'r')
	file1 = open("../gencode.v18.genes.patched_contigs.gtf_gene_tss", 'w')

	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file.readline()

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split()
		chr = line[0]
		tss = line[3]
		gene = line[9][1: -2]

		type = line[2]
		if type == 'transcript':
			file1.write(gene + '\t' + chr + '\t' + tss + '\n')

	file.close()
	file1.close()






	## get the list of all X, Y, MT genes
	file = open("../gencode.v18.genes.patched_contigs.gtf_gene_tss", 'r')
	file1 = open("../gencode.v18.genes.patched_contigs.gtf_gene_xymt", 'w')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		tss = line[2]

		if chr == 'X' or chr == 'Y' or chr == 'MT':
			file1.write(gene + '\n')
	file.close()
	file1.close()
