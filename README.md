# eQTL_script

All the data files should be in the upper folder of this directory, for appropriate processing.

The order of these scripts are:

** Genotype:

1. genotype\_ld\_prune.py

** Genotype backup (these should happen before above script):

1. script.py

2. post\_prune.py

3. [to-do] dosage\_qc\_chunk.py


** Expression:

eSample\_partition.py

gene\_preprocess.py

sample\_tissue\_preprocess.py

x. tissue\_hierarchy.py






Other scripts (not necessary in order):

1. gene\_tss\_prepare.py

2. GTEx\_beta\_extract.py




The pipeline for genotype QC and LD pruning:

As I may use different LD threshold (currently 0.5 for R^2) and the association threshold (currently 0.5 for R^2, the same with previous one) in the future, and there is hard drive usage issue, I record the procedure here.

(there are two scripts used here: “script.py" is the general script for processing all the chromosomes, following the below procedure; “post\_prune.py” is the one used for reversing the associations between pruned SNPs and the un-pruned SNPs; the two scripts are all in C2B2 “/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.5/44712/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v5.p1.c1.GRU/GenotypeFiles/phg000219.v4.GTEx\_Pilot\_Imputation.genotype-imputed-data.c1\_ld\_qc/")

tar zxvf “all.SNPs.tgz”
mkdir post\_prune
for all the chromosome “X”, do the following (step#04 — step#16):
tar zxvf “chrX.all.tgz”
cp chrX.tfam "chrX.tfam"
[QC] ./plink —tfile “chrX” —exclude “chrX.info4.maf05.exclusion.snplist.txt” —make-bed
[pruning] ./plink —bfile plink —indep-pairwise 50kb 5 0.5 (0.5 may possibly be adjusted into other values)
[LD statistics] ./plink --bfile plink --r2 --ld-snp-list plink.prune.out --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0.5 (0.5 can be adjusted into other values)
[LD statistics further] python post\_prune.py
mv plink.prune.in ./post\_prune/chr”X”.prune.in
mv plink.prune.out ./post\_prune/chr”X”.prune.out
mv plink.ld ./post\_prune/chr”X”.ld
mv chr”X".post\_prune.txt ./post\_prune/
rm chr22.*
rm data_imputed.*
rm plink.*

Then, I need to do the following procedure to get the genotype data (dosage) we can use in learning.

remove QC-ed SNPs from dosage file (simple; to be done with Python, in “/ifs/scratch/c2b2/ip_lab/sy2515/GTEx/data.v.5/44712/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v5.p1.c1.GRU/GenotypeFiles/phg000219.v4.GTEx_Pilot_Imputation.genotype-imputed-data.c1_dosage_qc”)
use the pruned SNP list to calculate the prior information for all un-pruned SNPs (if necessary), and get the potential for all un-pruned SNPs (will do with Python code; find the scheme first)
