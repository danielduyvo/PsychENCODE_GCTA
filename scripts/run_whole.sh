#!/bin/bash

# qsub options
#$ -w e
#$ -N ARG_NAME
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o ARG_NAME.log
#$ -e ARG_NAME.err
#$ -m a
#$ -M danieldu
#$ -t ARG_ARRAY

# Requires a complete GRM to be generated at data/${PROJECT}/output/grms/complete
# Can be generated with the following line using the file data/${PROJECT}/output/grm_chrs/all_chrs.txt
# gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/all_chrs.txt \
#     --out data/${PROJECT}/output/grms/complete

PROJECT="ARG_NAME"
THREADS=2
REDO=ARG_REDO

if $REDO; then
    GENE_ID=$(sed -n "${SGE_TASK_ID}s/\.hsq//gp" data/${PROJECT}/output/results/missing.txt)
    SGE_TASK_ID=$(awk "/${GENE_ID}/ {print NR}" data/${PROJECT}/input/phenotype_ids)
fi

LINE=$(sed -n ${SGE_TASK_ID}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR EXTRA
    printf "ID: %s\n" $ID

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -gt 0 ]]; then

            # Run GREML, defaulting to EM if AI fails
            gcta64 --reml --thread-num $THREADS --reml-alg 0 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --grm data/${PROJECT}/output/grms/complete\
                --pheno data/${PROJECT}/input/phenotype \
                --reml-lrt 1 \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --grm data/${PROJECT}/output/grms/complete\
                --pheno data/${PROJECT}/input/phenotype \
                --reml-lrt 1 \
                --out data/${PROJECT}/output/hsqs/$ID

        fi
    )
