# Heritability estimation with PsychENCODE data using GCTA

Here are scripts for running GCTA on the European frontal lobe samples from PsychENCODE. 
Tools and programs used inlude PLINK, GCTA and R.

## Folder structure
```
GCTA
|
+--data
|  |
|  \--PROJECT_NAME
|     |
|     +--output
|     |  |
|     |  +--fin_peds
|     |  |
|     |  +--graphs
|     |  |
|     |  +--grm_chrs
|     |  |
|     |  +--grm_ranges
|     |  |
|     |  +--grms
|     |  |
|     |  +--hsqs
|     |  |
|     |  +--mgrms
|     |  |
|     |  \--results
|     |
|     \--input
|        |
|        +--BED files
|        |
|        +--PED files
|        |
|        +--VCF files
|        |
|        +--Rdata files
|        |
|        +--Phenotype files for GREML
|        |
|        \--Removed related individuals
|
+--fin_genotype
|
+--fin_phenotype
|
+--scripts
|  |
|  +--plan.rmd
|  |
|  +--build_project.sh
|  |
|  +--run_bivariate_gcta.sh
|  |
|  \--run_gcta.sh for running GREML
|
\--README.md
```

## Usage

First time setting up, we'll need to generate necessary input files and set up the 
directory structure. The provided script `scripts/build_project.sh` will take in VCF 
and BED files and set up the initial directory hierarchy and necessary files for the 
GCTA analysis.

```
# Running analysis on European population, SNPs subsetted to those in HRC panel, gene 
# expression with sequencing PCs regressed out, calculating heritability from the SNPs 
# in a 1Mbase cis window:
# The following script will build a directory structure with the provided names, as 
# well as process and link the resulting files from the provided VCF and BED file.

./scripts/build_project.sh -f EUR gene_SPC HRC 1Mbase_cis genotype.vcf phenotype.bed
```

Afterwards, the provided script `scripts/run_gcta.sh` will run GCTA on the data.

```
# The -a flag allows us to choose the type of analysis
# The -t flag allows us to input a range of phenotypes we want to run GCTA with. 
# If we want to run against all provided phenotypes, we use 1-total number of phenotypes.
# The -w flag sets the cis window size.
# The -n flag points to the directory name that the previous script set up

./scripts/run_gcta.sh \
    -a cis \
    -n EUR_HRC_gene_SPC_1Mbase_cis \
    -t 1-$(cat fin_RNA_seq/EUR_gene_SPC/phenotype_ids | wc -l) \
    -w 1000000
```

Finally, we compile the HSQ output from GCTA into a single file, using 
`scripts/compile_results.sh`.

```
./scripts/compile_results.sh \
    -a cis \
    -n EUR_HRC_gene_SPC_1Mbase_cis
```

The results can be found at `data/EUR_HRC_gene_SPC_1Mbase_cis/output/results/all_variance.txt`
