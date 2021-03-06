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
                --reml-lrt 1 2 \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
                --remove data/${PROJECT}/input/removed_samples.txt \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --reml-lrt 1 2 \
                --out data/${PROJECT}/output/hsqs/$ID

            # Clean up files
            rm data/${PROJECT}/output/grm_chrs/$ID.txt
            rm data/${PROJECT}/output/mgrms/$ID.txt
            rm data/${PROJECT}/output/grms/${ID}_*

        fi
    )
