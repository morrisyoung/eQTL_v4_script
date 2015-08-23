## this is used after the sample processing script, and it will partition the eQTL samples into training set and testing set

import numpy as np

ratio = 0.75



if __name__ == '__main__':

	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples", 'r')
	file1 = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train", 'w')
	file2 = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test", 'w')
	
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		file1.write(tissue + '\t')
		file2.write(tissue + '\t')

		samples = line[1:]

		arr = np.arange(len(samples))
		np.random.shuffle(arr)

		i = 0
		while i < len(samples) * ratio:
			index = arr[i]
			sample = samples[index]
			file1.write(sample + '\t')
			i += 1
		file1.write('\n')

		while i < len(samples):
			index = arr[i]
			sample = samples[index]
			file2.write(sample + '\t')
			i += 1
		file2.write('\n')

	file.close()
	file1.close()
	file2.close()
