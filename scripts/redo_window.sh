#!/bin/bash

# qsub options
#$ -w e
#$ -N GCTA_analysis_window
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o window.log
#$ -e window.err
#$ -m a
#$ -M danieldu
#$ -t 1-142
PROJECT="EUR_SPC_HRC_1Mbase"
THREADS=2

SGE_TASK_ID=$((SGE_TASK_ID-1))
IDS=("ENSG00000066379" "ENSG00000132498" "ENSG00000154529" "ENSG00000155282" "ENSG00000156750" "ENSG00000170161" "ENSG00000170165" "ENSG00000170215" "ENSG00000170217" "ENSG00000170775" "ENSG00000172014" "ENSG00000182368" "ENSG00000183148" "ENSG00000184523" "ENSG00000185020" "ENSG00000185044" "ENSG00000187060" "ENSG00000188681" "ENSG00000196400" "ENSG00000196409" "ENSG00000196421" "ENSG00000196774" "ENSG00000197550" "ENSG00000198312" "ENSG00000198566" "ENSG00000204599" "ENSG00000204618" "ENSG00000204619" "ENSG00000204622" "ENSG00000204623" "ENSG00000204788" "ENSG00000204790" "ENSG00000204793" "ENSG00000204794" "ENSG00000204802" "ENSG00000204805" "ENSG00000204807" "ENSG00000204816" "ENSG00000204828" "ENSG00000204837" "ENSG00000206341" "ENSG00000206503" "ENSG00000212951" "ENSG00000212952" "ENSG00000215126" "ENSG00000215142" "ENSG00000215548" "ENSG00000216829" "ENSG00000219693" "ENSG00000223379" "ENSG00000223408" "ENSG00000223839" "ENSG00000224185" "ENSG00000224681" "ENSG00000224762" "ENSG00000225015" "ENSG00000225278" "ENSG00000225353" "ENSG00000225812" "ENSG00000225883" "ENSG00000226007" "ENSG00000226546" "ENSG00000227232" "ENSG00000227248" "ENSG00000227252" "ENSG00000227279" "ENSG00000227318" "ENSG00000227449" "ENSG00000227490" "ENSG00000227558" "ENSG00000227582" "ENSG00000227921" "ENSG00000228243" "ENSG00000229146" "ENSG00000229156" "ENSG00000229273" "ENSG00000229311" "ENSG00000229422" "ENSG00000230255" "ENSG00000230521" "ENSG00000230795" "ENSG00000230804" "ENSG00000230850" "ENSG00000230857" "ENSG00000230880" "ENSG00000231074" "ENSG00000231212" "ENSG00000231242" "ENSG00000231390" "ENSG00000231527" "ENSG00000231995" "ENSG00000232116" "ENSG00000232274" "ENSG00000232336" "ENSG00000232433" "ENSG00000232745" "ENSG00000232815" "ENSG00000232833" "ENSG00000232866" "ENSG00000233022" "ENSG00000233218" "ENSG00000233244" "ENSG00000233434" "ENSG00000233892" "ENSG00000233961" "ENSG00000234127" "ENSG00000234299" "ENSG00000234320" "ENSG00000234416" "ENSG00000234665" "ENSG00000234978" "ENSG00000235310" "ENSG00000235355" "ENSG00000235659" "ENSG00000235832" "ENSG00000236029" "ENSG00000236055" "ENSG00000236157" "ENSG00000236233" "ENSG00000236816" "ENSG00000237198" "ENSG00000237238" "ENSG00000237357" "ENSG00000237451" "ENSG00000237669" "ENSG00000237792" "ENSG00000237846" "ENSG00000237869" "ENSG00000238113" "ENSG00000238245" "ENSG00000238261" "ENSG00000240165" "ENSG00000240240" "ENSG00000240907" "ENSG00000241370" "ENSG00000242569" "ENSG00000242676" "ENSG00000243753" "ENSG00000255585" "ENSG00000260141" "ENSG00000261599" "ENSG00000263901" "ENSG00000264356" "ENSG00000264628" "ENSG00000268578" "ENSG00000269788" "ENSG00000270604" "ENSG00000270909" "ENSG00000271225")
SGE_TASK_ID=$(awk "/${IDS[$SGE_TASK_ID]}/ {print NR}" data/${PROJECT}/input/phenotype_ids)

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
            printf "data/%s/output/grms/%s_trans_int" $PROJECT $ID > data/${PROJECT}/output/grm_chrs/$ID.txt
            for ((i=1;i<=22;i++)); do
                if [[ $i -ne $CHR ]]; then
                    printf "\ndata/%s/output/grms/ped_file%s" $PROJECT $i >> data/${PROJECT}/output/grm_chrs/$ID.txt
                fi
            done

            # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
            printf "data/%s/output/grms/%s_cis\ndata/%s/output/grms/%s_trans" $PROJECT $ID $PROJECT $ID > \
                data/${PROJECT}/output/mgrms/$ID.txt

            # Split the cis and trans regions into separate PED files
            plink --bfile data/${PROJECT}/input/ped_file --extract range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_cis
            plink --bfile data/${PROJECT}/input/ped_file --exclude range data/${PROJECT}/output/grm_ranges/$ID.txt \
                --chr $CHR --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
                --out data/${PROJECT}/output/fin_peds/${ID}_trans_int

            # Generate the GRMs for the cis and trans regions
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_cis --make-grm-alg 0 \
                --out data/${PROJECT}/output/grms/${ID}_cis
            gcta64 --make-grm-bin --thread-num $THREADS --bfile data/${PROJECT}/output/fin_peds/${ID}_trans_int --make-grm-alg 0 --chr $CHR \
                --out data/${PROJECT}/output/grms/${ID}_trans_int
            gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/$ID.txt \
                --out data/${PROJECT}/output/grms/${ID}_trans

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
            rm data/${PROJECT}/output/grm_ranges/$ID.txt
            rm data/${PROJECT}/output/grm_chrs/$ID.txt
            rm data/${PROJECT}/output/mgrms/$ID.txt
            rm data/${PROJECT}/output/fin_peds/${ID}_*
            rm data/${PROJECT}/output/grms/${ID}_*

        fi
    )
