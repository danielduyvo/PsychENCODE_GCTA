#!/bin/bash

# qsub options
#$ -w e
#$ -N add_AFR_250kbase_GCTA_analysis_window
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o add_AFR_250kbase_window.log
#$ -e add_AFR_250kbase_window.err
#$ -m a
#$ -M danieldu
#$ -t 1-25774
SGE_TASK_ID=1000
PROJECT="AFR_SPC_HRC_250kbase"
# WINDOW=1000000 # 1Mbase
WINDOW=250000 # 250kbase
THREADS=2
LINE=$(sed -n ${SGE_TASK_ID}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR START END W_START W_END
    printf "ID: %s\tCHR: %s\n" $ID $CHR

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -gt 0 ]]; then

            # Now uses bash variable to decide window
            W_START=$(( ${START} - ${WINDOW} ))
            W_END=$(( ${END} + ${WINDOW} ))

            if [[ ${W_START} -lt 0 ]]; then
                W_START=0
            fi

            # Create the range file for PLINK to generate the cis region
            printf "%s %s %s R1" $CHR $W_START $W_END > data/${PROJECT}/output/grm_ranges/$ID.txt

            # Add number of cis SNPs to HSQ output
            plink --bfile data/${PROJECT}/input/ped_file \
                --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --remove data/${PROJECT}/input/removed_samples.txt --write-snplist \
                --out data/${PROJECT}/output/grm_ranges/$ID
            wc -l < data/${PROJECT}/output/grm_ranges/${ID}.snplist
            printf "cisSNPs\t%s" $(wc -l < data/${PROJECT}/output/grm_ranges/${ID}.snplist) \
                >> data/${PROJECT}/output/hsqs/$ID.hsq

            # Clean up files
            rm data/${PROJECT}/output/grm_ranges/$ID.*

        fi
    )
