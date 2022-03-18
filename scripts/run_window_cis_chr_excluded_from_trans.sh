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

PROJECT="ARG_NAME"
WINDOW=ARG_WINDOW
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

            # Now uses bash variable to decide window
            W_START=$(( ${START} - ${WINDOW} ))
            W_END=$(( ${END} + ${WINDOW} ))

            if [[ ${W_START} -lt 0 ]]; then
                W_START=0
            fi

            # Create the range file for PLINK to generate the cis and trans regions
            printf "%s %s %s R1" $CHR $W_START $W_END > data/${PROJECT}/output/grm_ranges/$ID.txt

            # Create file for GCTA pointing to which PED files belong to the trans region
            > data/${PROJECT}/output/grm_chrs/$ID.txt
            for ((i=1;i<=22;i++)); do
                if [[ $i -ne $CHR ]]; then
                    printf "\ndata/%s/output/grms/ped_file%s" $PROJECT $i >> data/${PROJECT}/output/grm_chrs/$ID.txt
                fi
            done

            # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
            printf "data/%s/output/grms/%s_cis\ndata/%s/output/grms/%s_trans" $PROJECT $ID $PROJECT $ID > \
                data/${PROJECT}/output/mgrms/$ID.txt

            # Split the cis into a separate PED file
            # If no cis SNPs, exit
            if plink --bfile data/${PROJECT}/input/ped_file --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_cis; then

                # Generate the GRMs for the cis and trans regions
                gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_cis --make-grm-alg 0 \
                    --out data/${PROJECT}/output/grms/${ID}_cis
                gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/$ID.txt \
                    --out data/${PROJECT}/output/grms/${ID}_trans

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

                # Add number of cis SNPs to HSQ output
                plink --bfile data/${PROJECT}/input/ped_file \
                    --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                    --remove data/${PROJECT}/input/removed_samples.txt --write-snplist \
                    --out data/${PROJECT}/output/grm_ranges/$ID
                wc -l < data/${PROJECT}/output/grm_ranges/${ID}.snplist
                printf "cisSNPs\t%s" $(wc -l < data/${PROJECT}/output/grm_ranges/${ID}.snplist) \
                    >> data/${PROJECT}/output/hsqs/$ID.hsq

                # Clean up files
                rm data/${PROJECT}/output/grms/${ID}_*
            fi
            # Clean up files
            rm data/${PROJECT}/output/grm_ranges/$ID.*
            rm data/${PROJECT}/output/grm_chrs/$ID.txt
            rm data/${PROJECT}/output/mgrms/$ID.txt
            rm data/${PROJECT}/output/fin_peds/${ID}_*

        fi
    )

