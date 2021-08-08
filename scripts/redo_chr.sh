#!/bin/bash

# qsub options
#$ -w e
#$ -N EUR_GCTA_analysis_chromosome
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o EUR_chr.log
#$ -e EUR_chr.err
#$ -m a
#$ -M danieldu
#$ -t 1-12

# This script splits the cis and trans region by chromosome instead
PROJECT="EUR_SPC_HRC_chrs"
THREADS=2
set -x
SGE_TASK_ID=$((SGE_TASK_ID-1))
IDS=("ENSG00000008128" "ENSG00000089693" "ENSG00000089723" "ENSG00000089737" "ENSG00000089775" "ENSG00000089902" "ENSG00000117151" "ENSG00000177225" "ENSG00000177234" "ENSG00000184178" "ENSG00000254131" "ENSG00000254141")
SGE_TASK_ID=$(awk "/${IDS[$SGE_TASK_ID]}/ {print NR}" data/${PROJECT}/input/phenotype_ids)

LINE=$(sed -n ${SGE_TASK_ID}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR START END W_START W_END
    printf "ID: %s\tCHR: %s\n" $ID $CHR

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -gt 0 ]]; then

            # Create file for GCTA pointing to which PED files belong to the trans region
            > data/${PROJECT}/output/grm_chrs/$ID.txt
            for ((i=1;i<=22;i++)); do
                if [[ $i -ne $CHR ]]; then
                    printf "\ndata/%s/output/grms/ped_file%s" $PROJECT $i >> data/${PROJECT}/output/grm_chrs/$ID.txt
                fi
            done

            # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
            printf "data/%s/output/grms/ped_file%s\ndata/%s/output/grms/%s_trans" $PROJECT $CHR $PROJECT $ID > \
                data/${PROJECT}/output/mgrms/$ID.txt

            # Generate the GRM for the trans region
            gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/$ID.txt \
                --out data/${PROJECT}/output/grms/${ID}_trans

            # Run GREML, defaulting to EM if AI fails
            gcta64 --reml --thread-num $THREADS --reml-alg 0 --reml-maxit 100 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID

            # Clean up files
            rm data/${PROJECT}/output/grm_chrs/$ID.txt
            rm data/${PROJECT}/output/mgrms/$ID.txt
            rm data/${PROJECT}/output/grms/${ID}_*

        fi
    )
