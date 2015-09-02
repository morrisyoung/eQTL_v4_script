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


	## get the x (185) individuals we will use in modeling
	file = open("../list_individual.txt", 'r')
	individual_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line
		individual_rep[individual] = 1
	file.close()


	## drop individual phenotypes that we will not use
	remove_list = []
	for individual in individual_var_rep:
		if individual not in individual_rep:
			remove_list.append(individual)
	for individual in remove_list:
		del individual_var_rep[individual]
	

	
	## to be used:
	##	individual_var_rep
	## to be generated:
	##	individual_batch_rep





	return individual_batch_rep






def batch_sample_extract():

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



	## get the y (185) individuals we will use in modeling
	file = open("../list_individual.txt", 'r')
	individual_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line
		individual_rep[individual] = 1
	file.close()





	return sample_batch_rep





if __name__ == '__main__':


	print "working on the batch..."

	batch_individual_extract()

	batch_sample_extract()

	print "done!"

