## this script is used to calculate the final prior score for all un-pruned SNPs, from prior scores from pruned and un-pruned SNPs and the association signals of pruned SNPs with un-pruend SNPs.


if __name__ == '__main__':

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









	##====== retriave the prior scores of all SNPs for these target tissues =======
	##====== get the association signals also =======
	##====== calculate the totalled score for all un-pruned SNPs =======
	##====== save them into file =======
	# under "../prior.score.unpruned/"











