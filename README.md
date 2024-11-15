# GWAS-Assisted and multi-trait genomic prediction for improvement of seed yield and canning quality traits in a black bean breeding panel

## Abstract
In recent years, black beans (Phaseolus vulgaris L.) have gained popularity in the U.S., with improved seed yield and canning quality being critical traits for new cultivars. Achieving genetic gains in these traits is often challenging due to negative trait associations and the need for specialized equipment and trained sensory panels for evaluation. This study investigates the integration of genomics and phenomics to enhance selection accuracy for these complex traits. We evaluated the prediction accuracy of single and multi-trait genomic prediction (GP) models, incorporating near-infrared spectroscopy (NIRS) data and markers identified through genome-wide association studies (GWAS). The models demonstrated moderate prediction accuracies for yield and canning appearance, and high accuracies for color retention. No significant differences were found between single-trait and multi-trait models within the same breeding cycle. However, across breeding cycles, multi-trait models outperformed single-trait models by up to 41% and 63% for canning appearance and seed yield, respectively. Interestingly, incorporating significant SNP markers identified by GWAS and NIRS data into the models tended to decrease prediction accuracy both within and between breeding cycles. As genotypes from the new breeding cycle were included, the models' prediction accuracy generally increased. Our findings underscore the potential of multi-trait models to enhance the prediction of complex traits such as seed yield and canning quality in dry beans and highlight the importance of continually updating the training dataset for effective GP implementation in dry bean breeding.

## Organization of this Repository

### Scripts:
0. R scripts for spatial analysis in the field to obtain BLUPs.
1. R scripts for comparing GBLUP and different kernels (K1, K2, K3, and KA).
2. R scripts for GWAS-assisted and multi-trait models across all traits.
3. R scripts for GWAS-assisted and multi-trait models specifically for canning appearance and yield.

### Data:
- GB_BLB.RData: An R list containing the phenotype, genotype, NIRs, and SNP positions.
- Raw Sequencing Data available at [http://www.ncbi.nlm.nih.gov/bioproject/1138671](https://www.ncbi.nlm.nih.gov/sra/PRJNA1138671)
