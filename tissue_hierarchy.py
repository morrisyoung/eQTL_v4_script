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

	print sample_tissue_map


	#file = open("", 'r')
