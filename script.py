## QC and prune SNPs from all chromosomes
import os


if __name__ == "__main__":

	## start here

	#chr_list = [22]
	chr_list = [22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

	for chr in chr_list:

		print "working on chr#",
		print chr

		##====== 1 ======
		print "extract the genotype data"
		os.system('tar zxvf chr' + str(chr) + '.all.tgz >/dev/null 2>&1')

		##====== 2 ======
		print "copy the chrX.tfam file"
		os.system('cp chrX.tfam chr' + str(chr) + '.tfam >/dev/null 2>&1')

		##====== 3 ======
		print "data quality control"
		os.system('./plink --tfile chr' + str(chr) + ' --exclude chr' + str(chr) + '.info4.maf05.exclusion.snplist.txt --make-bed >/dev/null 2>&1')

		##====== 4 ======
		print "LD pruning"
		os.system('./plink --bfile plink --indep-pairwise 50kb 5 0.5 >/dev/null 2>&1')

		##====== 5 ======
		print "LD statistics calculation"
		os.system('./plink --bfile plink --r2 --ld-snp-list plink.prune.out --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0.5 >/dev/null 2>&1')

		##====== 6 ======
		print "post-pruning processing"
		os.system('python post_prune.py >/dev/null 2>&1')


		##====== saving the files and cleaning ======
		print "saving results and cleaning"
		os.system('mv plink.prune.in ./post_prune/chr' + str(chr) + '.prune.in >/dev/null 2>&1')
		os.system('mv plink.prune.out ./post_prune/chr' + str(chr) + '.prune.out >/dev/null 2>&1')
		os.system('mv plink.ld ./post_prune/chr' + str(chr) + '.ld >/dev/null 2>&1')
		os.system('mv chr' + str(chr) + '.post_prune.txt ./post_prune/ >/dev/null 2>&1')
		os.system('rm chr' + str(chr) + '.* >/dev/null 2>&1')
		os.system('rm data_imputed.* >/dev/null 2>&1')
		os.system('rm plink.* >/dev/null 2>&1')




	print "done!"

