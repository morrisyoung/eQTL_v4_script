## extract and prepare the beta (empirical coefficients) from GTEx data analysis results
## TODO not finished here

if __name__ == '__main__':


	file_list = ["Adipose_Subcutaneous.portal.eqtl", "Muscle_Skeletal.portal.eqtl", "Artery_Aorta.portal.eqtl", "Nerve_Tibial.portal.eqtl", "Artery_Tibial.portal.eqtl", "Skin_Sun_Exposed_Lower_leg.portal.eqtl", "Esophagus_Mucosa.portal.eqtl", "Stomach.portal.eqtl", "Esophagus_Muscularis.portal.eqtl", "Thyroid.portal.eqtl", "Heart_Left_Ventricle.portal.eqtl", "Whole_Blood.portal.eqtl", "Lung.portal.eqtl"]

	for name in file_list:
		list = name.split('.')
		etissue = list[0]

		file = open("../GTEx_Analysis_V4_eQTLs/" + name, 'r')
		file1 = open("../GTEx_Analysis_V4_eQTLs/" + etissue + '.beta' , 'w')

		file.readline()

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')

			snp_id = line[0]
			snp_chr = line[1]
			snp_pos = line[2]
			gene_id = line[3]
			gene_pos = line[5]
			beta = line[7]
			print snp_id,
			print snp_chr,
			print snp_pos,
			print gene_id,
			print gene_pos,
			print beta
			print line
			break


		file.close()
		file1.close()
