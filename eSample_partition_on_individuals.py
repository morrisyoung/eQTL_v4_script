## this is used after the sample processing script, and it will partition the eQTL samples if individuals into training set and testing set (splitting individuals)

## this script will split the individuals into training set and testing set
## assuming each individual contains the same amount of samples from different tissues (on average), then all the samples are also 75%/25% split
## the output of this partition can be used in cross-tissue cis- multi-linear modeling



import numpy as np

ratio = 0.75



def sample_to_individual(sample):
	individual = sample[:9]
	return individual



if __name__ == '__main__':

	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples", 'r')

	tissue_sample_rep = {}		# element as list of samples (from that tissue)
	sample_tissue_map = {}		# element as tissue type (of that sample)
	individual_sample_rep = {}	# element as list of samples (from that individual)

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		sample_list = line[1:]

		## tissue_sample_rep
		tissue_sample_rep[tissue] = sample_list

		## sample_tissue_map
		for i in range(len(sample_list)):
			sample = sample_list[i]
			sample_tissue_map[sample] = tissue

		## individual_sample_rep
		for i in range(len(sample_list)):
			sample = sample_list[i]
			individual = sample_to_individual(sample)
			if individual not in individual_sample_rep:
				individual_sample_rep[individual] = [sample]
			else:
				individual_sample_rep[individual].append(sample)

	file.close()



	## pick up the training individuals and testing individuals, and save them
	file1 = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_individuals_train", 'w')
	file2 = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_individuals_test", 'w')

	individual_list = []
	for individual in individual_sample_rep:
		individual_list.append(individual)
	arr = np.arange(len(individual_list))
	np.random.shuffle(arr)


	i = 0
	while i < len(individual_list) * ratio:
		index = arr[i]
		individual = individual_list[index]
		file1.write(individual + '\t')
		## write the samples from this individual, to training set
		for j in range(len(individual_sample_rep[individual])):
			sample = individual_sample_rep[individual][j]
			file1.write(sample + '\t')
		file1.write('\n')
		i += 1


	while i < len(individual_list):
		index = arr[i]
		individual = individual_list[index]
		file2.write(individual + '\t')
		## write the samples from this individual, to testing set
		for j in range(len(individual_sample_rep[individual])):
			sample = individual_sample_rep[individual][j]
			file2.write(sample + '\t')
		file2.write('\n')
		i += 1

	file1.close()
	file2.close()

