## this script is used to extract the batch variables that can be utilized in the main routine
## 
##
##
## output:
##	batch_var_individual.txt
##	batch_var_sample.txt
##
##

import re




def DateToNum(date):
	## sample: "04/02/2012"
	date = date.split('/')
	month = int(date[0])
	day = int(date[1])
	year = int(date[2])

	## TODO: we start from 01/01/2010
	total = 0
	for i in range(2010, year):
		if i % 4 == 0:
			total += 366
		else:
			total += 365

	for i in range(1, month):
		if i in {1:1, 3:1, 5:1, 7:1, 8:1, 10:1, 12:1}:
			total += 31
		elif i == 2 and year % 4 == 0:
			total += 29
		elif i == 2 and year % 4 != 0:
			total += 28
		else:
			total += 30

	for i in range(1, day):
		total += 1

	return total





##================= issues for this routine ==================
## There are two types of batch variables: number-valued and string-valued; for number type, if we have any missing value, we won't use that variable, as there is no unbiased way to impute/quantify that missing number; for string type, we also quantify the missing value (they are already one special class in the original data), which simplifies the whole process. But not sure whether this is good enough. Under this condition, we remove 12 out of 172 batch variables.
def batch_individual_extract():

	individual_batch_var_list = []
	individual_batch_rep = {}

	## get the individual phenotype rep
	individual_batch_rep = {}
	file = open("../PhenotypeFiles/phs000424.v4.pht002742.v4.p1.c1.GTEx_Subject_Phenotypes.GRU.txt", 'r')
	individual_var_list = []
	individual_var_rep = {}
	count = 0
	while 1:
		line = file.readline()
		count += 1
		if count <= 10:
			continue
		if not line:
			break

		line = line[:-1]
		line = line.split('\t')
		if count == 11:
			individual_var_list = line[2:]
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
	


	##=========================== quantify all individual phenotypes (excluding some) ==============================
	## to be used:
	##	individual_var_list
	##	individual_var_rep
	## to be generated:
	##	individual_batch_var_list = []
	##	individual_batch_rep = {}
	##=========================== quantify all individual phenotypes (excluding some) ==============================
	## there are 'integer', 'decimal', 'enum_integer' and 'string' types
	## the 'integer' and 'decimal' types need to be linearly quantified, and 'enum_integer' and 'string' types can be directly quantified
	file = open("../PhenotypeFiles/phs000424.v4.pht002742.v4.p1.GTEx_Subject_Phenotypes.var_report.xml", 'r')
	count = 0
	var_type_rep = {}
	while 1:
		line = (file.readline()).strip()
		count += 1
		if not line:
			break
		## pattern: <variable id="phv00169061.v4.p1" var_name="SUBJID" calculated_type="string" reported_type="string">
		## may contain more than 1 item for each line
		it = re.finditer(r"<variable.*?>", line)
		for match in it:
			line1 = match.group()
			## get the var-type pair
			pattern = re.compile(r'var_name=\".*?\"')
			match = pattern.search(line1)
			var = match.group().split('=')[1][1:-1]
			pattern = re.compile(r'calculated_type=\".*?\"')
			match = pattern.search(line1)
			type = match.group().split('=')[1][1:-1]
			var_type_rep[var] = type
	file.close()

	##================== build the remove_list and remove_rep
	remove_list = ["TRCCLMP", "TRCHSTIN", "TRISCH", "DTHPRNINT"]  ## this is duplicated variables with others
	## extend the remove_list and remove_rep by checking missing batch value for all batch variables
	for index in range(len(individual_var_list)):
		var = individual_var_list[index]
		type = var_type_rep[var]
		count = 0
		for individual in individual_var_rep:
			value = individual_var_rep[individual][index]
			if value == '':
				count += 1
		if count > 0 and (type == 'integer' or type == 'decimal'):
			remove_list.append(var)
	remove_rep = {}
	for i in range(len(remove_list)):
		var = remove_list[i]
		remove_rep[var] = 1

	##==================== first of all, remove some individual phenotypes from 'remove_rep'
	remove_index = {}
	for i in range(len(individual_var_list)):
		var = individual_var_list[i]
		if var in remove_rep:
			remove_index[i] = var

	for i in range(len(individual_var_list)):
		var = individual_var_list[i]
		if i not in remove_index:
			individual_batch_var_list.append(var)
	for individual in individual_var_rep:
		individual_batch_rep[individual] = []
		for i in range(len(individual_var_rep[individual])):
			value = individual_var_rep[individual][i]
			if i not in remove_index:
				individual_batch_rep[individual].append(value)

	##==================== second, quantify all these batch variables
	individual_list = []
	for individual in individual_batch_rep:
		individual_list.append(individual)

	remove_list = []  ## if there are too many missing values, we should remove that batch variable
	for index in range(len(individual_batch_var_list)):
		var = individual_batch_var_list[index]
		type = var_type_rep[var]
		if type == 'integer' or type == 'decimal':
			value_list = []
			value_max = 0  ## TODO: the integer and decimal value will be scaled into minimum as 1
			for i in range(len(individual_list)):
				individual = individual_list[i]
				value = individual_batch_rep[individual][index]
				value = float(value)
				value_list.append(value)
				if value > value_max:
					value_max = value
			value_min = value_max
			for i in range(len(value_list)):
				value = value_list[i]
				if value < value_min:
					value_min = value
			for i in range(len(value_list)):
				individual = individual_list[i]
				individual_batch_rep[individual][index] = (value_list[i] - value_min) * 1.0 / value_max
		else:  ## 'enum_integer' or 'string' types: just quantify them
			value_list = []
			value_rep = {}
			count = 1
			for i in range(len(individual_list)):
				individual = individual_list[i]
				value = individual_batch_rep[individual][index]
				value_list.append(value)
				if value not in value_rep:
					value_rep[value] = count
					count += 1
			for value in value_rep:
				value_rep[value] = value_rep[value] * 1.0 / (count - 1)
			for i in range(len(value_list)):
				value = value_list[i]
				individual = individual_list[i]
				individual_batch_rep[individual][index] = value_rep[value]


	return (individual_batch_var_list, individual_batch_rep)






##================= issues for this routine ==================
## According to the same rule, we removed 3 out of 72 variables for sample batches.
def batch_sample_extract():

	sample_batch_var_list = []
	sample_batch_rep = {}

	## get the sample batch var
	sample_batch_rep = {}
	file = open("../PhenotypeFiles/phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt", 'r')
	sample_var_list = []
	sample_var_rep = {}
	count = 0
	while 1:
		line = file.readline()
		count += 1
		if count <= 10:
			continue
		if not line:
			break

		line = line[:-1]
		line = line.split('\t')
		if count == 11:
			sample_var_list = line[2:]
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


	##=========================== quantify all sample attributes (excluding some) ==============================
	## to be used:
	##	sample_var_list
	##	sample_var_rep
	## to be generated:
	##	sample_batch_var_list = []
	##	sample_batch_rep = {}
	##=========================== quantify all sample attributes (excluding some) ==============================
	## there are 'integer', 'decimal', 'enum_integer' and 'string' types
	## the 'integer' and 'decimal' types need to be linearly quantified, and 'enum_integer' and 'string' types can be directly quantified
	file = open("../PhenotypeFiles/phs000424.v4.pht002743.v4.p1.GTEx_Sample_Attributes.var_report.xml", 'r')
	var_type_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break
		## pattern: <variable id="phv00169061.v4.p1" var_name="SUBJID" calculated_type="string" reported_type="string">
		## may contain more than 1 item for each line
		it = re.finditer(r"<variable.*?>", line)
		for match in it:
			line1 = match.group()
			## get the var-type pair
			pattern = re.compile(r'var_name=\".*?\"')
			match = pattern.search(line1)
			var = match.group().split('=')[1][1:-1]
			pattern = re.compile(r'calculated_type=\".*?\"')
			match = pattern.search(line1)
			type = match.group().split('=')[1][1:-1]
			var_type_rep[var] = type
	file.close()

	##================== build the remove_list and remove_rep
	remove_list = ["SMGTC"]  ## this is duplicated variables with others
	## extend the remove_list and remove_rep by checking missing batch value for all batch variables
	for index in range(len(sample_var_list)):
		var = sample_var_list[index]
		type = var_type_rep[var]
		count = 0
		for sample in sample_var_rep:
			value = sample_var_rep[sample][index]
			if value == '':
				count += 1
		if count > 0 and (type == 'integer' or type == 'decimal'):
			remove_list.append(var)
	remove_rep = {}
	for i in range(len(remove_list)):
		var = remove_list[i]
		remove_rep[var] = 1

	##==================== first of all, remove some sample attributes from 'remove_rep'
	remove_index = {}
	for i in range(len(sample_var_list)):
		var = sample_var_list[i]
		if var in remove_rep:
			remove_index[i] = var

	for i in range(len(sample_var_list)):
		var = sample_var_list[i]
		if i not in remove_index:
			sample_batch_var_list.append(var)
	for sample in sample_var_rep:
		sample_batch_rep[sample] = []
		for i in range(len(sample_var_rep[sample])):
			value = sample_var_rep[sample][i]
			if i not in remove_index:
				sample_batch_rep[sample].append(value)

	##==================== second, quantify all these batch variables
	sample_list = []
	for sample in sample_batch_rep:
		sample_list.append(sample)

	remove_list = []  ## if there are too many missing values, we should remove that batch variable
	for index in range(len(sample_batch_var_list)):
		var = sample_batch_var_list[index]
		type = var_type_rep[var]
		if type == 'integer' or type == 'decimal':
			value_list = []
			value_max = 0  ## TODO: the integer and decimal value will be scaled into minimum as 1
			for i in range(len(sample_list)):
				sample = sample_list[i]
				value = sample_batch_rep[sample][index]
				value = float(value)
				value_list.append(value)
				if value > value_max:
					value_max = value
			value_min = value_max
			for i in range(len(value_list)):
				value = value_list[i]
				if value < value_min:
					value_min = value
			for i in range(len(value_list)):
				sample = sample_list[i]
				sample_batch_rep[sample][index] = (value_list[i] - value_min) * 1.0 / value_max
		else:  ## 'enum_integer' or 'string' types: just quantify them
			if var == "SMNABTCHD" or var == "SMGEBTCHD":
				value_list = []
				value_max = 0  ## TODO: the integer and decimal value will be scaled into minimum as 1
				for i in range(len(sample_list)):
					sample = sample_list[i]
					value = sample_batch_rep[sample][index]
					value = DateToNum(value)
					value_list.append(value)
					if value > value_max:
						value_max = value
				value_min = value_max
				for i in range(len(value_list)):
					value = value_list[i]
					if value < value_min:
						value_min = value
				for i in range(len(value_list)):
					sample = sample_list[i]
					sample_batch_rep[sample][index] = (value_list[i] - value_min) * 1.0 / value_max
			else:
				value_list = []
				value_rep = {}
				count = 1
				for i in range(len(sample_list)):
					sample = sample_list[i]
					value = sample_batch_rep[sample][index]
					value_list.append(value)
					if value not in value_rep:
						value_rep[value] = count
						count += 1
				for value in value_rep:
					value_rep[value] = value_rep[value] * 1.0 / (count - 1)
				for i in range(len(value_list)):
					value = value_list[i]
					sample = sample_list[i]
					sample_batch_rep[sample][index] = value_rep[value]


	return (sample_batch_var_list, sample_batch_rep)





if __name__ == '__main__':


	print "working on the batch extraction..."

	(individual_batch_var_list, individual_batch_rep) = batch_individual_extract()
	(sample_batch_var_list, sample_batch_rep) = batch_sample_extract()


	## save the results, for individual batches and sample batches
	## individual batch
	file = open("../batch_var_individual.txt", 'w')
	file.write('individual\t')
	for i in range(len(individual_batch_var_list)):
		var = individual_batch_var_list[i]
		file.write(var + '\t')
	file.write('\n')
	for individual in individual_batch_rep:
		file.write(individual + '\t')
		for i in range(len(individual_batch_rep[individual])):
			value = individual_batch_rep[individual][i]
			file.write(str(value) + '\t')
		file.write('\n')
	file.close()


	print len(sample_batch_var_list)


	## sample batch
	file = open("../batch_var_sample.txt", 'w')
	file.write('sample\t')
	for i in range(len(sample_batch_var_list)):
		var = sample_batch_var_list[i]
		file.write(var + '\t')
	file.write('\n')
	for sample in sample_batch_rep:
		file.write(sample + '\t')
		for i in range(len(sample_batch_rep[sample])):
			value = sample_batch_rep[sample][i]
			file.write(str(value) + '\t')
		file.write('\n')
	file.close()



	print "done!"
