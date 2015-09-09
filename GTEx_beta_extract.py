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

