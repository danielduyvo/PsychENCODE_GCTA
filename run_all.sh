#!/bin/bash

# qsub options
#$ -w e
#$ -N GCTA_analysis
#$ -l h_data=8G,h_rt=8:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o qsub.log
#$ -e qsub.err
#$ -m a
#$ -M danieldu
#$ -t 1-25774
PROJECT="1MBase"
THREADS=2
LINE=$(sed -n ${SGE_TASK_ID}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR START END W_START W_END
    printf "ID: %s\tCHR: %s\n" $ID $CHR

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -gt 0 ]]; then

            # Create the range file for PLINK to generate the cis and trans regions
            printf "%s %s %s R1" $CHR $W_START $W_END > data/${PROJECT}/output/grm_ranges/$ID.txt

            # Create file for GCTA pointing to which PED files belong to the trans region
            printf "grms/%s_trans_int" $ID > data/${PROJECT}/output/grm_chrs/$ID.txt
            for ((i=1;i<=22;i++)); do
                if [[ $i -ne $CHR ]]; then
                    printf "\ngrms/ped_file%s" $i >> data/${PROJECT}/output/grm_chrs/$ID.txt
                fi
            done

            # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
            printf "grms/%s_cis\ngrms/%s_trans" $ID $ID > data/${PROJECT}/output/mgrms/$ID.txt

            # Split the cis and trans regions into separate PED files
            plink --bfile ped_file --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_cis
            plink --bfile ped_file --exclude range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --chr $CHR --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_trans_int

            # Generate the GRMs for the cis and trans regions
            gcta64 --make-grm-bin --thread-num THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_cis --make-grm-alg 0 \
                --out data/${PROJECT}/output/grms/${ID}_cis
            gcta64 --make-grm-bin --thread-num THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_trans_int --make-grm-alg 0 --chr $CHR \
                --out data/${PROJECT}/output/grms/${ID}_trans_int
            gcta64 --make-grm-bin --thread-num THREADS --mgrm data/${PROJECT}/output/grm_chrs/$ID.txt \
                --out data/${PROJECT}/output/grms/${ID}_trans

            # Run GREML, defaulting to EM if AI fails
            gcta64 --reml --thread-num THREADS --reml-alg 0 --reml-maxit 100 --mpheno ${SGE_TASK_ID} \
                --mgrm data/${PROJECT}/output/mgrms/$ID.txt \
                --pheno data/${PROJECT}/input/phenotype \
                --out data/${PROJECT}/output/hsqs/$ID || \
            gcta64 --reml --thread-num THREADS --reml-alg 2 --reml-maxit 10000 --mpheno ${SGE_TASK_ID} \
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
