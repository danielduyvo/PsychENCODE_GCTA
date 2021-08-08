#!/bin/bash

# qsub options
#$ -w e
#$ -N EUR_GCTA_analysis_bksk
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o EUR_bksk.log
#$ -e EUR_bksk.err
#$ -m a
#$ -M danieldu
#$ -t 1-3

# Requires two additional GRMs:
# A complete GRM to be generated at data/${PROJECT}/output/grms/complete
# Can be generated with the following line using the file data/${PROJECT}/output/grm_chrs/all_chrs.txt
# gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/all_chrs.txt \
#     --out data/${PROJECT}/output/grms/complete
# A modified complete GRM with all unrelated individuals set to 0 at data/${PROJECT}/output/grms/bigk
# Can be generated with the following line using the file data/${PROJECT}/output/grms/complete
# gcta64 --grm data/${PROJECT}/output/grms/complete --thread-num $THREADS --make-bK 0.05 \
#     --out data/${PROJECT}/output/grms/bigk

PROJECT="EUR_SPC_HRC_bksk"
THREADS=2

SGE_TASK_ID=$((SGE_TASK_ID-1))
IDS=("ENSG00000249502" "ENSG00000249504" "ENSG00000249509")
SGE_TASK_ID=$(awk "/${IDS[$SGE_TASK_ID]}/ {print NR}" data/${PROJECT}/input/phenotype_ids)

LINE=$(sed -n ${SGE_TASK_ID}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR START END W_START W_END
    printf "ID: %s\n" $ID

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -gt 0 ]]; then

            # Create file for GCTA pointing to which GRM files to compare (Big K vs complete)
            printf "data/%s/output/grms/bigk\ndata/%s/output/grms/complete" $PROJECT $PROJECT > \
                data/${PROJECT}/output/mgrms/$ID.txt

            # Run GREML, defaulting to EM if AI fails
            gcta64 --reml --thread-num $THREADS --reml-alg 0 --reml-maxit 100 --mpheno ${SGE_TASK_ID} \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID

            # Clean up files
            rm data/${PROJECT}/output/mgrms/$ID.txt
        fi
    )