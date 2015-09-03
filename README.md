# eQTL_script

The C++ implementation for the learning program is at https://github.com/morrisyoung/eQTL_cplusplus.

All the data files should be in the upper folder of this directory, for appropriate processing.

## Genotype

1. genotype\_ld\_prune.py

### Genotype backup

(these should happen before above script, as this outputs what's used by above script; this is used in C2B2)

1. script\_qc\_ld\_pine.py
2. post\_prune.py
3. [to-do] dosage\_qc\_chunk.py

## Expression

1. sample\_tissue\_preprocess.py
2. gene\_preprocess.py
3. eSample\_partition.py
4. tissue\_hierarchy.py

## Other scripts

(for preprocessing information that may be used by the main routine, or mics)

1. gene\_tss\_prepare.py
2. GTEx\_beta\_extract.py
3. get\_individual.py
4. prior\_calculate.py
5. batch\_extract.py
6. practice\_cis\_detect.py
7. try.py


## The pipeline for genotype QC and LD pruning

As I may use different LD threshold (currently 0.5 for R^2) and the association threshold (currently 0.5 for R^2, the same with previous one) in the future, and there is hard drive usage issue, I record the procedure here.

(there are two scripts used here: “script.py" is the general script for processing all the chromosomes, following the below procedure; “post\_prune.py” is the one used for reversing the associations between pruned SNPs and the un-pruned SNPs; the two scripts are all in C2B2 “/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.5/44712/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v5.p1.c1.GRU/GenotypeFiles/phg000219.v4.GTEx\_Pilot\_Imputation.genotype-imputed-data.c1\_ld\_qc/")

```
1. tar zxvf “all.SNPs.tgz”
2. mkdir post\_prune
3. for all the chromosome “X”, do the following (step#04 — step#16):
4. tar zxvf “chrX.all.tgz”
5. cp chrX.tfam "chrX.tfam"
6. [QC] ./plink —tfile “chrX” —exclude “chrX.info4.maf05.exclusion.snplist.txt” —make-bed
7. [pruning] ./plink —bfile plink —indep-pairwise 50kb 5 0.5 (0.5 may possibly be adjusted into other values)
8. [LD statistics] ./plink --bfile plink --r2 --ld-snp-list plink.prune.out --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0.5 (0.5 can be adjusted into other values)
9. [LD statistics further] python post\_prune.py
10. mv plink.prune.in ./post\_prune/chr”X”.prune.in
11. mv plink.prune.out ./post\_prune/chr”X”.prune.out
12. mv plink.ld ./post\_prune/chr”X”.ld
13. mv chr”X".post\_prune.txt ./post\_prune/
14. rm chr22.*
15. rm data_imputed.*
16. rm plink.\*
```

Then, I need to do the following procedure to get the genotype data (dosage) we can use in learning.

1. Remove QC-ed SNPs from dosage file (simple; to be done with Python, in “/ifs/scratch/c2b2/ip\_lab/sy2515/GTEx/data.v.5/44712/PhenoGenotypeFiles/RootStudyConsentSet\_phs000424.GTEx.v5.p1.c1.GRU/GenotypeFiles/phg000219.v4.GTEx\_Pilot\_Imputation.genotype-imputed-data.c1\_dosage\_qc”).


## Expression data processing

1. eTissue is defined as GTEx tissues that have >= 60 effective samples (having genotype information).
2. non-Null gene is defined as "at least \portion of the eSamples have rpkm value >= \threshold", where \portion is 0.5 and \threshold is 0.1 currently.
3. We randomly select 75% of all eSamples in each eTissue as the training set, and the left as the testing set. The training set and testing set are prepared before learning and testing in the main routine.
4. We have the hierarchichal clustering results for fully processed expression file (sample dimension, gene dimension), but we have two versions, one for normalized expression matrix (quantile), and another for un-normalized.


## Fold enrichment of chromatin states

(learned from Roadmap Epigenomics project) on GWAS SNPs (data downloaded from ENCODE project)

We start from the first paper of this series: “http://www.nature.com/encode/threads/impact-of-functional-information-on-understanding-variation".

From that paper, we use the following file as the GWAS SNPs (originally from NHGRI GWAS SNP catalog June 2011): Gwascatalog.june2011.positions.bed, from “http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/supplementary/integration_data_jan2011/byDataType/GWAS/jan2011/".

I treat each snp site as one source, though there may be several snp-phenotype associations at that site. This assumes that the associated phenotypes may be correlated themselves, so that site only contributes once. That’s why we use the above file, "Gwascatalog.june2011.positions.bed".

As we also know the learned chromatin states with their chromosome positions (from Roadmap Epigenomics project), we can calculate the following: the % of the GWAS SNP set overlapping with one chromatin state divided by the % of the total segments this chromatin state class makes up. This is the enrichment value we need for down-stream analysis.

I remove all information (learned chromatin states and GWAS SNPs) from X and Y chromatin, as currently we are only interested in the autosome cis- prior for genotypes, and we actually don't have SNPs from X and Y in GTEx project.

There are some ready-to-use results from Roadmap people, in 2010, "Discovery and characterization of chromatin states for systematic annotation of the human genome", in "http://www.nature.com/nbt/journal/v28/n8/fig_tab/nbt.1662_F4.html". But as there are only 1,640 GWAS SNPs (other than 4492 we are using here) considered, we will not use this smaller version of analysis.

Here is a mapping from Roadmap Epigenomics project epigenome to the tissue type we have in GTEx project (reversed) if there is:

```
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
```

## Batch extraction

There are two types of batch variables: number-valued and string-valued; for number type, if we have any missing value, we won’t use that variable, as there is no unbiased way to impute/quantify that missing number; for string type, we also quantify the missing value (they are already one special class in the original data), which simplifies the whole process. But not sure whether this is good enough. Under this condition, we remove 12 out of 172 individual batch variables, and 3 out of 72 sample batch variables.

The tables we refered are as followed:

#### Individual phenotypes
ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v4.p1/pheno_variable_summaries/phs000424.v4.pht002742.v4.GTEx_Subject_Phenotypes.data_dict.xml
ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v4.p1/pheno_variable_summaries/phs000424.v4.pht002742.v4.p1.GTEx_Subject_Phenotypes.var_report.xml

#### sample attributes
ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v4.p1/pheno_variable_summaries/phs000424.v4.pht002743.v4.GTEx_Sample_Attributes.data_dict.xml
ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v4.p1/pheno_variable_summaries/phs000424.v4.pht002743.v4.p1.GTEx_Sample_Attributes.var_report.xml


## Statistics and what we have

### Statistics

* Currently we have **185 individuals**, **824,176 SNPs** (dosage data), **17 eTissues**, **1,889 eSamples**, **20,603 non-Null genes**. We can load them all into memory.

### Genotype relevant

* We know the lists of **pruned SNPs** and **un-pruned SNPs** (due to LD pruning, after QC), in "../genotype\_185\_dosage\_matrix\_qc/post\_prune/chrX\_prune.in" or "../genotype\_185\_dosage\_matrix\_qc/post\_prune/chrX\_prune.out".
* We have the **dosage data** for all un-pruned SNPs in "../genotype\_185\_dosage\_matrix\_qc/chrX/...", and we also have their name and position in the same folder.
* We have the **beta** (only significant association) from GTEx project, which we can utilize in the initialization of our learning, in "../GTEx\_Analysis\_V4\_eQTLs/...".
* We know the **association coefficients (R^2)** of pruned SNPs with their representative un-pruned SNPs, in "../genotype\_185\_dosage\_matrix\_qc/post\_prune/chrX.post\_prune.txt".
* We have the **enrichment value of chromatin states for all pruned and un-pruned SNPs**, in "../prior.score/etissue#/...". Some eTissues don't have this information, as the GTEx tissues are not fully consistent with the Roadmap Epigenomics tissues. We use "../prior.tissue.epigenome.map" to **map the eTissues in GTEx to epigenomics in Roadmap**, and "../prior.tissue.index.map" to map eTissues to an index for convenience of saving the enrichment values.


### Expression relevant

* We have **fully processed expression data** (remove samples having no genotype and not in eTissues; remove Null genes; quantile normalize across all genes), as "../GTEx\_Data\_2014-01-17\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_rpkm.gct\_processed\_2\_gene\_normalized".
* We have the **mean expression level** (after fully processing, normalized or not) as "../GTEx\_Data\_2014-01-17\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_rpkm.gct\_processed\_2\_gene\_normalized\_tissue\_mean", or "../GTEx\_Data\_2014-01-17\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_rpkm.gct\_processed\_2\_gene\_tissue\_mean". We also have the learned **tissue hierarchy** as "../tissue\_hierarchy\_normalized.txt" or "../tissue\_hierarchy\_unnormalized.txt" correspondingly.
* We have **gene TSS** as "../gencode.v18.genes.patched\_contigs.gtf\_gene\_tss".
* As X, Y or MT genes don't have cis- SNPs (we only have autosome genotypes from GTEx), but we still consider them in our framework (they may contribute to some cell env variables), we have **list of X, Y and MT genes**.
* We have the **eSample list of each eTissue** (sample size >= 60), as "../phs000424.v4.pht002743.v4.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type\_60\_samples", and we further partition them into **training set** as "../phs000424.v4.pht002743.v4.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type\_60\_samples\_train" and **testing set** as "phs000424.v4.pht002743.v4.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type\_60\_samples\_test".


