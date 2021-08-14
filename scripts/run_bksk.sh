#!/bin/bash

# qsub options
#$ -w e
#$ -N GCTA_analysis_bksk
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o bksk.log
#$ -e bksk.err
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
# gcta64 --grm data/${PROJECT}/output/grms/complete --thread-num $THREADS --make-bK 0.025 \
#     --out data/${PROJECT}/output/grms/bigk

PROJECT="EUR_SPC_HRC_bksk_0.025"
THREADS=2

REDO=0

if [[ $REDO -eq 0 ]]; then
    GENE_ID=$(sed -n "${SGE_TASK_ID}s/\.hsq//gp" data/${PROJECT}/output/results/missing.txt)
    SGE_TASK_ID=$(awk "/${GENE_ID}/ {print NR}" data/${PROJECT}/input/phenotype_ids)
fi

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
                --reml-lrt 1 2 \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --reml-lrt 1 2 \
                --out data/${PROJECT}/output/hsqs/$ID

            # Clean up files
            rm data/${PROJECT}/output/mgrms/$ID.txt
        fi
    )
