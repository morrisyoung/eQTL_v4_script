## this script is used to extract the batch variables that can be utilized in the main routine
## 
##
##
## output:
##	batch_var_individual.txt
##	batch_var_sample.txt
##
##


def batch_individual_extract():

	individual_batch_var_list = []
	individual_batch_rep = {}

	## get the individual phenotype rep
	individual_batch_rep = {}
	file = open("../PhenotypeFiles/phs000424.v4.pht002742.v4.p1.c1.GTEx_Subject_Phenotypes.GRU.txt", 'r')
	var_list = []
	individual_var_rep = {}
	count = 0
	while 1:
		line = file.readline()
		count += 1
		if count <= 10:
			continue
		if not line:
			break


		line = line.split('\t')
		if count == 11:
			var_list = line[2:]
		else:
			individual = line[1]
			individual_var_rep[individual] = line[2:]
	file.close()


	## get the x (currently 185) individuals we will use in modeling
	file = open("../list_individual.txt", 'r')
	individual_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line
		individual_rep[individual] = 1
	file.close()


	## drop individuals that we will not use
	remove_list = []
	for individual in individual_var_rep:
		if individual not in individual_rep:
			remove_list.append(individual)
	for individual in remove_list:
		del individual_var_rep[individual]
	

	
	## to be used:
	##	individual_var_rep
	## to be generated:
	##	individual_batch_var_list = []
	##	individual_batch_rep = {}








	return (individual_batch_var_list, individual_batch_rep)






def batch_sample_extract():

	sampme_batch_var_list = []
	sample_batch_rep = {}

	## get the sample batch var
	sample_batch_rep = {}
	file = open("../PhenotypeFiles/phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt", 'r')
	var_list = []
	sample_var_rep = {}
	count = 0
	while 1:
		line = file.readline()
		count += 1
		if count <= 10:
			continue
		if not line:
			break

		line = line.split('\t')
		if count == 11:
			var_list = line[2:]
		else:
			sample = line[1]
			sample_var_rep[sample] = line[2:]
	file.close()


	## get the y (currently 1889) eSamples we will use in modeling
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples", 'r')
	sample_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = (line.split('\t'))[1:]
		for sample in line:
			sample_rep[sample] = 1
	file.close()


	## drop samples that we will not use
	remove_list = []
	for sample in sample_var_rep:
		if sample not in sample_rep:
			remove_list.append(sample)
	for sample in remove_list:
		del sample_var_rep[sample]




	## to be used:
	##	sample_var_rep
	## to be generated:
	##	sampme_batch_var_list = []
	##	sample_batch_rep = {}










	return (sampme_batch_var_list, sample_batch_rep)





if __name__ == '__main__':


	print "working on the batch extraction..."

	(individual_batch_var_list, individual_batch_rep) = batch_individual_extract()
	(sampme_batch_var_list, sample_batch_rep) = batch_sample_extract()

	print "done!"

