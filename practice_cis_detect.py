## for one specific gene, detect the association of all cis- SNPs, in blood sample.

import time
#from scipy import stats
import numpy as np



if __name__ == '__main__':


	time_start = time.time()

	####=== build the sample-tissue repository ===####
	file = open("./data/sample_tissue.data", "r")
	sample_tissue = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		if len(line) != 2: ## there is missing data
			continue
		sample = line[0]
		tissue = line[1]
		sample_tissue[sample] = tissue
	file.close()

	print "done sample_tissue..."
	####==========================================####


	####=== get the 185 individuals that will be used in eQTL analysis ===####
	file = open("./data/individuals_185.data", "r")
	individuals = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		individual = line[0] ## like "GTEX-TKQ2"; this should be called subjectID
		index = int(line[1]) ## one number within 185
		individuals[individual] = index
	file.close()

	print "done individuals..."
	####==========================================####



	####=== get all the dosage information (genotypes) for each SNP for all individuals on Chr22 ===####
	## we only get the SNPs for 185 individuals on Chr22
	file = open("./data/dosage_185.data", "r")
	dosage = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		SNP = line[0]
		pos = int(line[1])
		dosage[SNP] = [pos]
		dosage[SNP].extend(map(lambda x: float(x), line[2:]))
	file.close()

	print "done dosage..."
	####==========================================####



	####=== get the position of each gene ===####
	## get all the genes (transcripts) with their positions only on Chr22
	file = open("./data/gene_pos_Chr22.data", "r")
	gene_pos = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		gene = line[0]
		start = int(line[1])
		end = int(line[2])
		gene_pos[gene] = (start, end)
	file.close()

	print "done gene_pos..."
	####==========================================####



	####=== get the expression files for all samples (from different individuals accross tissues) ===####
	## we only detect the association for genes in Chr22, so here we only get the expression for all samples for genes only on Chr22
	file = open("./data/gene_expression.data", "r")
	sample_list = []
	gene_expression = {}
	count = 0
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		if count == 0: ## this is the sample list
			for sample in line[1:]:
				sample_list.append(sample)
			count += 1
			continue

		gene = line[0]
		gene_expression[gene] = map(lambda x: float(x), line[1:])
		count += 1
	file.close()

	print "done gene_expression..."
	####==========================================####


	time_end = time.time()
	print "finish preparing the raw data..."
	print "time used is",
	print time_end - time_start,
	print "seconds..."
	print "Now start working!!..."
	## now we have sample_tissue{}, individuals{}, dosage{}, gene_pos{}, sample_list[] and gene_expression{}



	####==========================####
	####=== start working here ===####
	####==========================####
	for gene in gene_expression:  # or select by hand
		print "the gene we are interested in is:",
		print gene
		

		##==== now get the expression file for this specific gene ====
		# sample_list
		expression_array = gene_expression[gene]


		##==== now extract all the samples from blood tissue, from expression_array[] with sample_list ====
		sample_blood = []
		expression_blood = []
		for i in range(len(sample_list)):
			sample = sample_list[i]
			if sample_tissue[sample] == "Blood" and sample[:9] in individuals: ## this sample is from blood, and we have the genotype information
				sample_blood.append(sample)
				expression_blood.append(expression_array[i])
		print "we get",
		print len(sample_blood),
		print "samples from blood tissue, for which we have genotype information."



		##==== now try to get all the cis- SNPs for this specific gene ====
		gene_start = gene_pos[gene][0]
		gene_end = gene_pos[gene][1]
		SNP_cis = {}
		for SNP in dosage:
			pos = dosage[SNP][0]
			if pos - gene_start < 1000000 and pos - gene_start > -1000000: ## this is a cis- SNP
				SNP_cis[SNP] = dosage[SNP]
		print "we have",
		print len(SNP_cis),
		print "that we are interested in as cis- SNPs."


		##==== now get the genotypes for all this individuals of the samples ====
		dosage_blood_index = []  # find the index first
		for sample in sample_blood:
			subject = sample[:9]
			index = individuals[subject]
			dosage_blood_index.append(index)
		genotype_blood = {}
		for SNP in SNP_cis:
			genotype_blood[SNP] = []
			for index in dosage_blood_index:
				genotype_blood[SNP].append(SNP_cis[SNP][index])

		###=================================####
		### now we have the following:
		##	--- gene
		##	--- expression_blood
		##	--- genotype_blood
		### need to find out the associations
		###=================================####


		p_list = []  # store all the cis- SNPs and their P values
		for SNP in genotype_blood:
			genotype = genotype_blood[SNP]
			expression = expression_blood

			## least-squares regression			
			x = np.array(genotype)
			y = np.array(expression)
			slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

			p_list.append((SNP, p_value))


		print "We can get the following significant associated SNPs (with their P values):"
		p_threshold = 0.01
		for pair in p_list:
			if pair[1] < p_threshold:
				print pair[0],
				print pair[1]



		break ## only detect one single gene's expression here
