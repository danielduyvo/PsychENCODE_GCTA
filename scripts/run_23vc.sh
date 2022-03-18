#!/bin/bash

# qsub options
#$ -w e
#$ -N ARG_NAME
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o ARG_NAME.log
#$ -e ARG_NAME.err
#$ -m a
#$ -M danieldu
#$ -t ARG_ARRAY


# This script splits the cis and trans region by chromosome instead

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
    read ID CHR START END EXTRA
    printf "ID: %s\tCHR: %s\n" $ID $CHR


    # Create file for GCTA pointing to which PED files belong to the trans region
    > data/${PROJECT}/output/mgrms/$ID.txt
    for ((i=1;i<=22;i++)); do
        printf "data/%s/output/grms/ped_file%s\n" $PROJECT $i >> \
            data/${PROJECT}/output/mgrms/$ID.txt
    done

    # Run GREML, defaulting to EM if AI fails
    gcta64 --reml --thread-num $THREADS --reml-alg 0 --reml-maxit 100 --mpheno ${SGE_TASK_ID} \
        --remove data/${PROJECT}/input/removed_samples.txt \
        --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
        --pheno data/${PROJECT}/input/phenotype \
        --reml-lrt 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 \
        --out data/${PROJECT}/output/hsqs/$ID || \
    gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
        --remove data/${PROJECT}/input/removed_samples.txt \
        --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
        --pheno data/${PROJECT}/input/phenotype \
        --reml-lrt 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 \
        --out data/${PROJECT}/output/hsqs/$ID

    # Clean up files
    rm data/${PROJECT}/output/mgrms/$ID.txt

    )

