## this is used to calculate the prior score (chromatin state enrichment score for GWAS SNPs)
## the input files are as followed:
##	"../prior.all.mnemonics.bedFiles/...", which contains all the chromatin states for each epigenome (each tissue type)
##	"../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_60", which contains the eTissue list
##	"../prior.tissue.epigenome.map", which contains the mapping from eTissues to the epigenomes (if there is)
##	"../genotype_185_dosage_matrix_qc/post_prune/chr#.prune.in|out", which contains all the analyzed genotypes (pruned or un-pruned)
##	"../genotype_185_dosage_matrix/chr1/SNP_info.txt", which contains the snp positions of all snps (before QC)
## the output of this script:
##	prior scores grouped by eTissues and chromosomes


## more:
##	1. we should mark all eTissues with numbers, or "etissue#", for convenience of saving all the learned score
##	2. we should save the learned scores grouped by tissue types and chromosomes: "../prior.score/etissue#i/chr#j.in.score" and "../prior.score/etissue#i/chr#j.out.score"
##	3. in all the analysis, we should remove all the X, Y chromosomes, because we don't have SNPs (GTEx) from these chromosomes

import os




def snp_cs_map(chr_rep):
	## map the snps (pruned and un-pruned in GTEx after QC) into their chromatin states
	## this indeed gives us the enrichment value of each snp, and we should save that into files

	## we can use binary search to get the position (chromatin state) of one snp

	for i in range(22):
		# what we need:
		prune_in_list = []
		prune_out_list = []
		snp_info_rep = {}
		chr = "chr" + str(i+1)
		
		# fill in them
		file = open("../genotype_185_dosage_matrix_qc/post_prune/" + chr + ".prune.in", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			prune_in_list.append(line)
		file.close()
		file = open("../genotype_185_dosage_matrix_qc/post_prune/" + chr + ".prune.out", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			prune_out_list.append(line)
		file.close()

		file = open("../genotype_185_dosage_matrix/" + chr + "/SNP_info.txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')
			snp = line[0]
			pos = int(line[1])
			snp_info_rep[snp] = pos
		file.close()

		
		##===== now calculate =====
		# fill in them:
		state_in_list = []
		state_out_list = []
		for i in range(len(prune_in_list)):
			snp = prune_in_list[i]
			pos = snp_info_rep[snp]

			index = len(chr_rep[chr][0])/2  # the current pointer
			amount = len(chr_rep[chr][0])/4  # the amount to be shifted
			state = ''
			while 1:
				if index == (len(chr_rep[chr][0])-1):
					state = chr_rep[chr][1][-1]
					break
				if index == 0:
					state = chr_rep[chr][1][0]
					break
				if chr_rep[chr][0][index] <= pos and chr_rep[chr][0][index+1] > pos:
					state = chr_rep[chr][1][index]
					break
				elif chr_rep[chr][0][index] > pos:
					index = index - amount
					amount = amount / 2
					if amount == 0:
						amount = 1
				else:
					index = index + amount
					amount = amount / 2
					if amount == 0:
						amount = 1
			state_in_list.append(state)

		for i in range(len(prune_out_list)):
			snp = prune_out_list[i]
			pos = snp_info_rep[snp]

			index = len(chr_rep[chr][0])/2  # the current pointer
			amount = len(chr_rep[chr][0])/4  # the amount to be shifted
			state = ''
			while 1:
				if index == (len(chr_rep[chr][0])-1):
					state = chr_rep[chr][1][-1]
					break
				if index == 0:
					state = chr_rep[chr][1][0]
					break
				if chr_rep[chr][0][index] <= pos and chr_rep[chr][0][index+1] > pos:
					state = chr_rep[chr][1][index]
					break
				elif chr_rep[chr][0][index] > pos:
					index = index - amount
					amount = amount / 2
					if amount == 0:
						amount = 1
				else:
					index = index + amount
					amount = amount / 2
					if amount == 0:
						amount = 1
			state_out_list.append(state)



		print len(state_in_list)
		print len(state_out_list)

		## transform into prior value
		prior_in_list = []
		prior_out_list = []
		







	return





def cs_portion_gwas_portion(index):
	## given the epigenome and the GWAS SNPs, calculate:
	## the portion of chromosome of each chromatin state, and the portion of overlapping GWAS SNPs of each chromatin states
	## the ratio is indeed the enrichment value of that chromatin states

	##================== the portion of chromosome of each state ========================
	filename = "../prior.all.mnemonics.bedFiles/" + index + "_15_coreMarks_mnemonics.bed"
	file = open(filename, 'r')

	state_rep = {}  ## TODO: this is one of what we need finally
	length_total = 0

	chr_rep = {}  ## key as the chr#, and value is a list of two lists: start_list, state_list

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		chr = line[0]
		start = int(line[1])
		state = line[3]

		if chr == 'chrX' or chr == 'chrY':
			continue
		
		if chr in chr_rep:
			chr_rep[chr][0].append(start)
			chr_rep[chr][1].append(state)
		else:
			chr_rep[chr] = [[start], [state]]


		dist = int(line[2]) - int(line[1])
		if state in state_rep:
			state_rep[state] += dist
			length_total += dist
		else:
			state_rep[state] = dist
			length_total += dist
	file.close()

	for state in state_rep:
		state_rep[state] = state_rep[state] * 1.0 / length_total

	''' test
	print "there are",
	print len(state_rep),
	print "states overall."
	total = 0
	for state in state_rep:
		print state,
		print state_rep[state]
		total += state_rep[state]
	print "they are summed to",
	print total
	'''


	##================== the portion of overlapping GWAS SNPs of each chromatin states
	file = open("../Gwascatalog.june2011.positions.bed", 'r')
	gwas_rep = {}  # key as the chromatin state, and value as the number of gwas snps
	num_total = 0
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		chr = line[0]
		start = int(line[1])

		if chr == 'chrX' or chr == 'chrY':
			continue
	
		num_total += 1

		index = len(chr_rep[chr][0])/2  # the current pointer
		amount = len(chr_rep[chr][0])/4  # the amount to be shifted
		state = ''
		while 1:
			if index == (len(chr_rep[chr][0])-1):
				state = chr_rep[chr][1][-1]
				break
			if index == 0:
				state = chr_rep[chr][1][0]
				break
			if chr_rep[chr][0][index] <= start and chr_rep[chr][0][index+1] > start:
				state = chr_rep[chr][1][index]
				break
			elif chr_rep[chr][0][index] > start:
				index = index - amount
				amount = amount / 2
				if amount == 0:
					amount = 1
			else:
				index = index + amount
				amount = amount / 2
				if amount == 0:
					amount = 1

		if state in gwas_rep:
			gwas_rep[state] += 1
		else:
			gwas_rep[state] = 1

	file.close()
	for state in gwas_rep:
		gwas_rep[state] = gwas_rep[state] * 1.0 / num_total


	# test
	'''
	print "there are",
	print len(gwas_rep),
	print "states overall."
	total = 0
	for state in gwas_rep:
		print state,
		print gwas_rep[state]
		total += gwas_rep[state]
	print "they are summed to",
	print total
	'''


	##=============== the enrichment value ================
	# state_rep, gwas_rep
	enrich_rep = {}
	for state in state_rep:
		enrich_rep[state] = gwas_rep[state] / state_rep[state]


	# test
	'''
	print "there are",
	print len(enrich_rep),
	print "states overall."
	total = 0
	for state in enrich_rep:
		print state,
		print enrich_rep[state]
		total += enrich_rep[state]
	print "they are summed to",
	print total
	'''


	

	snp_cs_map(chr_rep, enrich_rep)






	return








if __name__ == "__main__":

	

	print "test starts..."


	#for i in range(1, 18):
	#	os.system("mkdir ../prior.score/etissue" + str(i))


	index = "E042"
	cs_portion_gwas_portion(index)


	print "test done..."


