## for doing the ld pruning for qc-ed data

## the input: 1. chunked dosage data (by chromosome and individual); 2. the pruned.in list from the PLINK pruning program
## the output: chunked dosage data (by chromosome and individual) with only un-pruned data in

## TODO for new dataset: the individual ID may contain more than 9 characters




if __name__ == '__main__':



	##========= get 185 individuals first ===========
	file = open("../GTEx_5M_185_Dec2012.Eigenvectors.txt", 'r')
	file.readline()
	rep_individual = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line[0:9]
		rep_individual[individual] = 1
	file.close()



	for chr in range(1, 23):



		##================ get all SNPs after QC =================
		rep_snp = {}

		file = open("../genotype_185_dosage_matrix_qc/post_prune/chr" + str(chr) + ".prune.in", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			snp = line
			rep_snp[snp] = 1
		file.close()

		print "working on chr#",
		print chr,
		print "and there are",
		print len(rep_snp),
		print "un-pruned SNPs."


		##========== get the positions of all SNPs left after QC ==========
		file = open("../genotype_185_dosage_matrix_qc/chr" + str(chr) + "/SNP_info.txt", 'r')
		rep_index = {}  # record the positions of all SNPs after QC
		index = 0
		while 1:
			line = (file.readline()).strip()
			index += 1
			if not line:
				break

			line = line.split(" ")
			snp = line[0]
			if snp in rep_snp:
				rep_index[index] = 1

		file.close()


		##========== re-write all dosage file ============
		for individual in rep_individual:
			filename = "../genotype_185_dosage_matrix_qc/chr" + str(chr) + "/SNP_dosage_" + individual + ".txt"

			## go through all the snps for this individual
			file = open(filename, 'r')
			dosage_list = []
			index = 0
			while 1:
				line = (file.readline()).strip()
				index += 1
				if not line:
					break

				dosage = line
				if index in rep_index:
					dosage_list.append(dosage)
			file.close()

			## write the list to the file
			file = open(filename, 'w')
			for dosage in dosage_list:
				file.write(dosage + '\n')
			file.close()



		##========== re-write all SNP_info.txt files ============
		filename = "../genotype_185_dosage_matrix_qc/chr" + str(chr) + "/SNP_info.txt"

		## go through all the snps for this individual
		file = open(filename, 'r')
		info_list = []
		index = 0
		while 1:
			line = (file.readline()).strip()
			index += 1
			if not line:
				break

			info = line
			if index in rep_index:
				info_list.append(info)
		file.close()

		## write the list to the file
		file = open(filename, 'w')
		for info in info_list:
			file.write(info + '\n')
		file.close()
