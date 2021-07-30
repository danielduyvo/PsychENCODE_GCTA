# Heritability estimation with PsychENCODE data using GCTA

Here are scripts for running GCTA on the European frontal lobe samples from PsychENCODE. 
Tools and programs used inlude PLINK, GCTA and R.

## File structure
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
+--scripts
|  |
|  +--plan.rmd
|  |
|  \--run_all.sh
|
\--README.md
```

