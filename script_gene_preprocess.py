## this can be used either before or after the script_sample_tissue_preprocess.py
## this is used to remove NULL genes (according to the )

if __name__ == '__main__':

	# start from here
	file = oepn("./GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file.readline()
	file.readline()
	file.readline()

	while 1:
		break


	file.close()

	
