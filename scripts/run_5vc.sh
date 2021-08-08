#!/bin/bash

# qsub options
#$ -w e
#$ -N AFR_SPC_HRC_bksk_1Mbase
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o AFR_SPC_HRC_bksk_1Mbase.log
#$ -e AFR_SPC_HRC_bksk_1Mbase.err
#$ -m a
#$ -M danieldu
#$ -t 1-25774

PROJECT="AFR_SPC_HRC_bksk_1Mbase"
WINDOW=1000000 # 1Mbase
# WINDOW=250000 # 250kbase
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

            # Create the range file for PLINK to generate the cis and trans regions
            printf "%s %s %s R1" $CHR $W_START $W_END > data/${PROJECT}/output/grm_ranges/$ID.txt

            # Create file for GCTA pointing to which PED files belong to the trans region
            printf "data/%s/output/grms/%s_trans_int" $PROJECT $ID > data/${PROJECT}/output/grm_chrs/$ID.txt
            for ((i=1;i<=22;i++)); do
                if [[ $i -ne $CHR ]]; then
                    printf "\ndata/%s/output/grms/ped_file%s" $PROJECT $i >> data/${PROJECT}/output/grm_chrs/$ID.txt
                fi
            done

            # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
            printf "data/%s/output/grms/%s_cis_bigk\n" $PROJECT $ID > \
                data/${PROJECT}/output/mgrms/$ID.txt
            printf "data/%s/output/grms/%s_cis\n" $PROJECT $ID >> \
                data/${PROJECT}/output/mgrms/$ID.txt
            printf "data/%s/output/grms/%s_trans_bigk\n" $PROJECT $ID >> \
                data/${PROJECT}/output/mgrms/$ID.txt
            printf "data/%s/output/grms/%s_trans\n" $PROJECT $ID >> \
                data/${PROJECT}/output/mgrms/$ID.txt

            # Split the cis and trans regions into separate PED files
            plink --bfile data/${PROJECT}/input/ped_file --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --make-bed --out data/${PROJECT}/output/fin_peds/${ID}_cis

            plink --bfile data/${PROJECT}/input/ped_file --exclude range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --chr $CHR --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_trans_int

            # Generate the GRMs for the cis and trans regions, split between big K small K
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_cis \
                --make-bK 0.025 --make-grm-alg 0 \
                --out data/${PROJECT}/output/grms/${ID}_cis_bigk
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_trans_int \
                --make-bK 0.025 --make-grm-alg 0 --chr $CHR \
                --out data/${PROJECT}/output/grms/${ID}_trans_int_bigk
            gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/$ID.txt \
                --make-bk 0.025 --out data/${PROJECT}/output/grms/${ID}_trans_bigk
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_cis \
                --make-grm-alg 0 --out data/${PROJECT}/output/grms/${ID}_cis
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_trans_int \
                --make-grm-alg 0 --chr $CHR \
                --out data/${PROJECT}/output/grms/${ID}_trans_int
            gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/$ID.txt \
                --out data/${PROJECT}/output/grms/${ID}_trans

            # Run GREML, defaulting to EM if AI fails
            gcta64 --reml --thread-num $THREADS --reml-alg 0 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID

            # Clean up files
            rm data/${PROJECT}/output/grm_ranges/$ID.txt
            rm data/${PROJECT}/output/grm_chrs/$ID.txt
            rm data/${PROJECT}/output/mgrms/$ID.txt
            rm data/${PROJECT}/output/fin_peds/${ID}_*
            rm data/${PROJECT}/output/grms/${ID}_*

        fi
    )
