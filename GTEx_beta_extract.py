## extract and prepare the beta (empirical coefficients) from GTEx data analysis results

## what we will do in this script:
### extract the beta from source files
### map the old file name to the tissue name (its index) used in main program
### save the beta to the new file name
### save all the tissues (with their indices) in a separate list, which we can refer easily later on

### I will use a tissue index for mapping all GTEx tissues to a number; that is only for convenience, and it's not related to tissue indexing in other folders

### later on, if there are new significant associations from different eTissues, we only need to add them into filename_map


size_cis_window = 1000000


if __name__ == '__main__':


	##======== get the tissue mapping, and the tissue-index mapping ==========
	filename_map = {"Adipose_Subcutaneous.portal.eqtl": "Adipose - Subcutaneous",
			"Muscle_Skeletal.portal.eqtl": "Muscle - Skeletal",
			"Artery_Aorta.portal.eqtl": "Artery - Aorta",
			"Nerve_Tibial.portal.eqtl": "Nerve - Tibial",
			"Artery_Tibial.portal.eqtl": "Artery - Tibial",
			"Skin_Sun_Exposed_Lower_leg.portal.eqtl": "Skin - Sun Exposed (Lower leg)",
			"Esophagus_Mucosa.portal.eqtl": "Esophagus - Mucosa",
			"Stomach.portal.eqtl": "Stomach",
			"Esophagus_Muscularis.portal.eqtl": "Esophagus - Muscularis",
			"Thyroid.portal.eqtl": "Thyroid",
			"Heart_Left_Ventricle.portal.eqtl": "Heart - Left Ventricle",
			"Whole_Blood.portal.eqtl": "Whole Blood",
			"Lung.portal.eqtl": "Lung"}

	tissue_index_map = {}
	index = 1
	file = open("../GTEx_Analysis_V4_eQTLs/etissue_list.txt", 'w')
	for filename in filename_map:
		tissue = filename_map[filename]
		tissue_index_map[tissue] = index
		file.write(tissue + "\t" + str(index) + '\n')
		index += 1
	file.close()




	##======== get the pruning association map ========
	rep_pruned_unpruned = {}
	rep_all = {}
	for chr in range(1, 23):
		file = open("../genotype_185_dosage_matrix_qc/post_prune/chr" + str(chr) + ".post_prune.txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split("\t")
			unpruned = line[0]
			rep_all[unpruned] = 1
			for i in range(1, len(line)):
				pair = line[i].split(" ")
				pruned = pair[0]
				r2 = float(pair[1])
				rep_pruned_unpruned[pruned] = (unpruned, r2)
				rep_all[pruned] = 1


	print "*"

	##======== extracting the beta from original files, and saving them into indexed filenames ==========
	for filename in filename_map:
		tissue = filename_map[filename]
		index = tissue_index_map[tissue]
		filename_new = "etissue" + str(index) + ".beta"

		file1 = open("../GTEx_Analysis_V4_eQTLs/" + filename, 'r')
		file2 = open("../GTEx_Analysis_V4_eQTLs/" + filename_new, 'w')
		file1.readline()  # remove the ID line

		while 1:
			line = (file1.readline()).strip()
			if not line:
				break

			line = line.split('\t')

			snp_id = line[0]
			#snp_chr = line[1]
			snp_pos = line[2]
			gene_id = line[3]
			gene_pos = line[5]
			beta = line[7]
			'''
			print snp_id,
			print snp_chr,
			print snp_pos,
			print gene_id,
			print gene_pos,
			print beta
			print line
			'''
			
			if snp_id not in rep_all:
				continue

			if snp_id in rep_pruned_unpruned:
				unpruned = rep_pruned_unpruned[snp_id][0]
				r2 = rep_pruned_unpruned[snp_id][1]
				snp_id = unpruned
				beta = str(float(beta) * r2)

			# test the distance between SNP and gene (whether cis-? should be)
			pos_snp = long(snp_pos)
			pos_chr = long(gene_pos)
			if abs(pos_snp - pos_chr) > size_cis_window:
				print snp_id,
				print gene_id,
				print beta
				continue

			file2.write(gene_id + "\t" + snp_id + "\t" + beta + "\n")

		file1.close()
		file2.close()

