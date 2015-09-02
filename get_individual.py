## this is used to get the individual list from some source file.


if __name__ == "__main__":

	file = open("../GTEx_5M_185_Dec2012.Eigenvectors.txt", 'r')
	file1 = open("../list_individual.txt", 'w')
	file.readline()

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		individual = line[0:9]
		file1.write(individual + '\n')

	file.close()
	file1.close()
