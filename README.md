# eQTL_script

All the data files should be in the upper folder of this directory, for appropriate processing.

The order of these scripts are:

==============================================

** Genotype:

1. genotype\_ld\_prune.py

** Genotype backup (these should happen before above script):

1. script.py

2. post\_prune.py

3. [to-do] dosage\_qc\_chunk.py

==============================================

** Expression:

1. sample\_tissue\_preprocess.py

2. eSample\_partition.py

3. gene\_preprocess.py

4. tissue\_hierarchy.py

==============================================

Other scripts (not necessary in order):

1. gene\_tss\_prepare.py

2. GTEx\_beta\_extract.py

3. practice\_cis\_detect.py



==============================================
==============================================
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
rm plink.\*

Then, I need to do the following procedure to get the genotype data (dosage) we can use in learning.

remove QC-ed SNPs from dosage file (simple; to be done with Python, in “/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.5/44712/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v5.p1.c1.GRU/GenotypeFiles/phg000219.v4.GTEx\_Pilot\_Imputation.genotype-imputed-data.c1\_dosage\_qc”)
use the pruned SNP list to calculate the prior information for all un-pruned SNPs (if necessary), and get the potential for all un-pruned SNPs (will do with Python code; find the scheme first)


==============================================
==============================================
Fold enrichment of chromatin states (learned from Roadmap Epigenomics project) on GWAS SNPs (data downloaded from ENCODE project):

We start from the first paper of this series: “http://www.nature.com/encode/threads/impact-of-functional-information-on-understanding-variation".

From that paper, we use the following file as the GWAS SNPs (originally from NHGRI GWAS SNP catalog June 2011): Gwascatalog.june2011.positions.bed, from “http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/supplementary/integration\_data\_jan2011/byDataType/GWAS/jan2011/".

I treat each snp site as one source, though there may be several snp-phenotype associations at that site. This assumes that the associated phenotypes may be correlated themselves, so that site only contributes once. That’s why we use the above file, "Gwascatalog.june2011.positions.bed".

As we also know the learned chromatin states with their chromosome positions (from Roadmap Epigenomics project), we can calculate the following: the % of the GWAS SNP set overlapping with one chromatin state divided by the % of the total segments this chromatin state class makes up. This is the enrichment value we need for down-stream analysis.

There are some ready-to-use results from Roadmap people, in 2010, "Discovery and characterization of chromatin states for systematic annotation of the human genome", in "http://www.nature.com/nbt/journal/v28/n8/fig\_tab/nbt.1662\_F4.html". But as there are only 1,640 GWAS SNPs considered, we will not use this smaller version of analysis.

Here is a mapping from Roadmap Epigenomics project epigenome to the tissue type we have in GTEx project (reversed) if there is:


Whole Blood:
E062
E034
E045
E033
E044
E043
E039
E041
E042
E040
E037
E048
E038
E047
E029
E031
E035
E051
E050
E036
E032
E046
E030
E115
E116
E123
E124

Cells - Transformed fibroblasts:
E126
E128

Thyroid:
None

Esophagus - Mucosa:
E079

Pancreas:
E087
E098

Artery - Tibial:
None

Testis:
None

Adipose - Subcutaneous:
E025
E063

Nerve - Tibial:
None

Artery - Aorta:
E065

Stomach:
E111
E092
E110
E094

Colon - Transverse:
E076
E106
E075

Esophagus - Muscularis:
E079

Skin - Sun Exposed (Lower leg):
E055
E056
E059
E061
E057
E058
E126
E127

Heart - Left Ventricle:
E095

Muscle - Skeletal:
E108
E107
E120
E121

Lung:
E096



==============================================
==============================================
Stats:

1. currently we have 824176 SNPs used in the modeling.

2. We know the lists of pruned SNPs and un-pruned SNPs (due to LD pruning, after QC), in "../genotype\_185\_dosage\_matrix\_qc/post\_prune/chrX\_prune.in" or "../genotype\_185\_dosage\_matrix\_qc/post\_prune/chrX\_prune.out".

3. We know the association coefficients (R^2) of pruned SNPs with their representative un-pruned SNPs, in "../genotype\_185\_dosage\_matrix\_qc/post\_prune/chrX.post\_prune.txt".

4. We have the dosage data for all un-pruned SNPs in "../genotype\_185\_dosage\_matrix\_qc/chrX/...", and we also have their name and position in the same folder.

5. We have the beta (only significant association) from GTEx project, which we can utilize in the initialization of our learning, in "../GTEx\_Analysis\_V4\_eQTLs/...".

6. We have fully processed expression data (remove samples having no genotype and not in eTissues; remove Null genes; quantile normalize across all genes), as "../GTEx\_Data\_2014-01-17\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_rpkm.gct\_processed\_2\_gene\_normalized".

7. We have the mean expression level (after fully processing, normalized or not) as "../GTEx\_Data\_2014-01-17\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_rpkm.gct\_processed\_2\_gene\_normalized\_tissue\_mean", or "../GTEx\_Data\_2014-01-17\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_rpkm.gct\_processed\_2\_gene\_tissue\_mean". We also have the learned tissue hierarchy as "../tissue\_hierarchy\_normalized.txt" or "../tissue\_hierarchy\_unnormalized.txt" correspondingly.

8. We have gene TSS as "../gencode.v18.genes.patched\_contigs.gtf\_gene\_tss".

9. As X, Y or MT genes don't have cis- SNPs (we only have autosome genotypes from GTEx), but we still consider them in our framework (they may contribute to some cell env variables), we have list of X, Y and MT genes.

10. We have the sample list of each eTissues (sample size >= 60), as "../phs000424.v4.pht002743.v4.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type", and we further partition them into training set as "../phs000424.v4.pht002743.v4.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type\_60\_samples\_train" and testing set as "phs000424.v4.pht002743.v4.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type\_60\_samples\_test".

