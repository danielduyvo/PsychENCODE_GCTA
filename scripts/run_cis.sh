#!/bin/bash

# qsub options
#$ -w e
#$ -N GCTA_analysis_cis
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o cis.log
#$ -e cis.err
#$ -m a
#$ -M danieldu
#$ -t 1-25774
PROJECT="EUR_SPC_HRC_cis"
THREADS=2
LINE=$(sed -n ${SGE_TASK_ID}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR START END W_START W_END
    printf "ID: %s\tCHR: %s\n" $ID $CHR

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -gt 0 ]]; then

            # Create the range file for PLINK to generate the cis region
            printf "%s %s %s R1" $CHR $W_START $W_END > data/${PROJECT}/output/grm_ranges/$ID.txt

            # Split the cis and trans regions into separate PED files
            plink --bfile data/${PROJECT}/input/ped_file --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_cis

            # Generate the GRMs for the cis region
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_cis --make-grm-alg 0 \
                --out data/${PROJECT}/output/grms/${ID}_cis

            # Run GREML, defaulting to EM if AI fails
            gcta64 --reml --thread-num $THREADS --reml-alg 0 --reml-maxit 100 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --grm data/${PROJECT}/output/grms/${ID}_cis\
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --grm data/${PROJECT}/output/grms/${ID}_cis\
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID

            # Clean up files
            rm data/${PROJECT}/output/fin_peds/${ID}_*
            rm data/${PROJECT}/output/grms/${ID}_*

        fi
    )
