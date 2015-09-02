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


def snp_cs_map(chr_rep, enrich_rep):
	## map the snps (pruned and un-pruned in GTEx after QC) into their chromatin states, and then enrichment value
	list1 = []
	list2 = []

	for i in range(22):
		# what we need:
		prune_in_list = []
		prune_out_list = []
		snp_info_rep = {}
		chr = "chr" + str(i+1)

		#print "chr#",
		#print chr
		
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
		# record the states of all snps in these two lists
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


		'''
		print "in and out have",
		print len(state_in_list),
		print "and",
		print len(state_out_list),
		print "snps on this chromosome"
		'''

		## transform into prior value, and prepare to return these results
		prior_in_list = []
		prior_out_list = []
		for state in state_in_list:
			enrich = enrich_rep[state]
			prior_in_list.append(enrich)
		for state in state_out_list:
			enrich = enrich_rep[state]
			prior_out_list.append(enrich)

		
		list1.append(prior_in_list)
		list2.append(prior_out_list)


	return (list1, list2)



def cs_portion_gwas_portion(index):
	## given the epigenome and the GWAS SNPs, calculate:
	## the portion of chromosome of each chromatin state, and the portion of overlapping GWAS SNPs of each chromatin states
	## the ratio is indeed the enrichment value of that chromatin states
	##================== the portion of chromosome of each state ========================
	filename = "../prior.all.mnemonics.bedFiles/" + index + "_15_coreMarks_mnemonics.bed"
	file = open(filename, 'r')

	state_rep = {}  ## key as the state, and value as its portion of the total length of chromosome
	length_total = 0

	chr_rep = {}  ## key as the chr#, and value is a list of two lists: start_list, state_list
			## this is to be used later on for binary finding of the chromatin state of one SNP

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


	print "got chr portion of each state"

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

	print "got gwas snp portion of each state"


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
		if state not in gwas_rep:
			enrich_rep[state] = 0
		else:
			enrich_rep[state] = gwas_rep[state] / state_rep[state]

	print "got enrichment value for each state"

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

	print "working on all GTEx snps (pruned or un-pruned)..."

	(list1, list2) = snp_cs_map(chr_rep, enrich_rep)
	return (list1, list2)





if __name__ == "__main__":

	print "start..."
	'''
	for i in range(1, 18):
		os.system("mkdir ../prior.score/etissue" + str(i))
	'''
	tissue_rep = {}
	file = open("../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_60", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		tissue = line.split('\t')[0]
		tissue_rep[tissue] = []
	file.close()

	file = open("../prior.tissue.epigenome.map", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		if len(line) == 1:
			continue

		tissue = line[0]
		epi_list = line[1].split(' ')
		tissue_rep[tissue] = epi_list
	file.close()

	print "we have the eTissues and Epigenomes mapping as following:",
	print tissue_rep

	## build the "prior.tissue.index.map"
	file = open("../prior.tissue.index.map", 'w')
	tissue_index_rep = {}
	index = 0
	for tissue in tissue_rep:
		index = index + 1
		file.write(tissue + '\t' + str(index) + '\n')
		tissue_index_rep[tissue] = index
	file.close()

	for tissue in tissue_rep:
		if len(tissue_rep[tissue]) == 0:
			print "no epigenomes for GTEx tissue",
			print tissue
			continue
	
		print "working on GTEx tissue",
		print tissue

		index = tissue_index_rep[tissue]

		L1 = []
		L2 = []

		for epigenome in tissue_rep[tissue]:
			(list1, list2) = cs_portion_gwas_portion(epigenome)  # list1 is for prune.in and list2 is for prune.out
			if len(L1) == 0:
				L1 = list1[:]
				L2 = list2[:]
			else:
				for i in range(22):
					for j in range(len(L1[i])):
						L1[i][j] = L1[i][j] + list1[i][j]
					for j in range(len(L2[i])):
						L2[i][j] = L2[i][j] + list2[i][j]
		for i in range(22):
			for j in range(len(L1[i])):
				L1[i][j] = L1[i][j] / len(tissue_rep[tissue])
			for j in range(len(L2[i])):
				L2[i][j] = L2[i][j] / len(tissue_rep[tissue])

		for chr in range(1, 23):
			file = open("../prior.score/etissue" + str(index) + "/chr" + str(chr) + ".in.score", 'w')
			for i in range(len(L1[chr-1])):
				enrich = L1[chr-1][i]
				file.write(str(enrich) + '\n')
			file.close()

			file = open("../prior.score/etissue" + str(index) + "/chr" + str(chr) + ".out.score", 'w')
			for i in range(len(L2[chr-1])):
				enrich = L2[chr-1][i]
				file.write(str(enrich) + '\n')
			file.close()

	print "done!"


