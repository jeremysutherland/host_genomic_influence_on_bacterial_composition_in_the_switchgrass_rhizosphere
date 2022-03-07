# host_genomic_influence_on_bacterial_composition_in_the_switchgrass_rhizosphere

This repository contains 9 scripts used in the journal article "Host Genomic Influence on Bacterial Composition in the Switchgrass Rhizosphere"

map.R contains the R code used to generate the map for Figure 1

phyloseq.R contains the R code for all of the phyloseq analyses including code for Figure 2, Supplemental Figure 2, Supplemental Figure 3, Supplemental Figure 4, and Supplemental Table 1.

corrs.R contains the R code used to generate the scatterplots for Figure 3

h2.R contains the R code used to calculate narrow-sense heritability and generate Figure 4

RDA.R contains the R code for the RDA analysis and Figure 5.

statgenGWAS.R contains the R code used for the genome wide association study and to generate Figure 6 and Supplemental Figure 5. 

fix_snipe.py contains the python code used to modify the HapMap v2 SNP matrix (https://doi.org/10.5061/dryad.mp6cp)

sommer_MAF_snp_mat.R contains the R code used to further process the HapMap v2 SNP matrix and filter SNPs.

dada2.R contains the R code for processing the raw reads deposited in the NCBI SRA under BioProject PRJNA689762
