## this script is used to calculate the final prior score for all un-pruned SNPs, from prior scores from pruned and un-pruned SNPs and the association signals of pruned SNPs with un-pruend SNPs.

import os


if __name__ == '__main__':

	'''
	## etissueX
	for i in range(1, 23):
		print i
		os.system("mkdir ../prior.score.unpruned/etissue" + str(i))
	'''

	##====== first get which tissue has the prior knowledge =======
	# from "../prior.tissue.epigenome.map"
	file = open("../prior.tissue.epigenome.map", 'r')
	prior_tissue_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		if len(line) == 1:
			continue

		tissue = line[0]
		prior_tissue_rep[tissue] = 1
	file.close()

	file = open("../prior.tissue.index.map", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		tissue = line[0]
		index = line[1]
		
		if tissue in prior_tissue_rep:
			prior_tissue_rep[tissue] = index
	file.close()

	##====== and save them into file =======
	# under "../prior.score.unpruned/"

	file = open("../prior.score.unpruned/prior_tissue_index.txt", 'w')
	for tissue in prior_tissue_rep:
		index = prior_tissue_rep[tissue]
		file.write(tissue + '\t' + index + '\n')
	file.close()

	##====== retriave the prior scores of all SNPs for these target tissues =======
	##====== get the association signals also =======
	##====== calculate the totalled score for all un-pruned SNPs =======
	##====== save them into file =======
	# under "../prior.score.unpruned/"
	# with prior_tissue_rep
	# with "../prior.score/etissueX/chrY.in.score" or "../prior.score/etissueX/chrY.out.score", and "../genotype_185_dosage_matrix_qc/post_prune/chrY.post_prune.txt"

	for tissue in prior_tissue_rep:
		index = prior_tissue_rep[tissue]
		print "we have tissue",
		print tissue,
		print "which has index as",
		print index

		for chr in range(1, 23):
			snp_in_list = []
			snp_out_list = []

			prior_in_list = []
			prior_out_list = []

			prior_in_rep = {}
			prior_out_rep = {}

			file = open("../genotype_185_dosage_matrix_qc/post_prune/chr" + str(chr) + ".prune.in", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				snp = line
				snp_in_list.append(snp)
			file.close()

			file = open("../genotype_185_dosage_matrix_qc/post_prune/chr" + str(chr) + ".prune.out", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				snp = line
				snp_out_list.append(snp)
			file.close()

			file = open("../prior.score/etissue" + index + "/chr" + str(chr) + ".in.score", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				prior = line
				prior_in_list.append(float(prior))
			file.close()

			file = open("../prior.score/etissue" + index + "/chr" + str(chr) + ".out.score", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				prior = line
				prior_out_list.append(float(prior))
			file.close()

			for i in range(len(snp_in_list)):
				snp = snp_in_list[i]
				prior = prior_in_list[i]
				prior_in_rep[snp] = prior

			for i in range(len(snp_out_list)):
				snp = snp_out_list[i]
				prior = prior_out_list[i]
				prior_out_rep[snp] = prior

			# what we may use: snp_in_list, prior_out_rep, prior_in_rep

			assoc_rep = {}
			file = open("../genotype_185_dosage_matrix_qc/post_prune/chr" + str(chr) + ".post_prune.txt", 'r')
			while 1:
				line = (file.readline()).strip()
				if not line:
					break

				line = line.split('\t')
				unpruned = line[0]
				assoc_rep[unpruned] = []
				for i in range(1, len(line)):
					pair = line[i].split(' ')
					pruned = pair[0]
					r2 = float(pair[1])
					assoc_rep[unpruned].append((pruned, r2))

			# what we may use: snp_in_list, prior_out_rep, prior_in_rep, assoc_rep

			prior_final_list = []
			for i in range(len(snp_in_list)):
				snp = snp_in_list[i]
				prior_final = prior_in_rep[snp]
				if snp in assoc_rep:
					for j in range(len(assoc_rep[snp])):
						pair = assoc_rep[snp][j]
						pruned = pair[0]
						r2 = pair[1]
						prior = prior_out_rep[pruned]
						prior_final += r2 * prior
				prior_final_list.append(prior_final)



			# save the results, for this eTissue, on this chromosome

			file = open("../prior.score.unpruned/etissue" + str(index) + "/chr" + str(chr) + ".score", 'w')
			for i in range(len(snp_in_list)):
				prior_final = prior_final_list[i]
				file.write(str(prior_final) + '\n')
			file.close()
