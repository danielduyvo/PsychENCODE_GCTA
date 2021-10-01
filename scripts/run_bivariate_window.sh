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

PROJECT="ARG_NAME"
WINDOW=ARG_WINDOW

GENE="ARG_GROUP"
TRANSCRIPTS=$(wc -l < data/${PROJECT}/input/${GENE}_ids)
THREADS=2

mkdir -p data/$PROJECT/output/fin_peds/${GENE}
mkdir -p data/$PROJECT/output/graphs/${GENE}
mkdir -p data/$PROJECT/output/grm_chrs/${GENE}
mkdir -p data/$PROJECT/output/grm_ranges/${GENE}
mkdir -p data/$PROJECT/output/grms/${GENE}
mkdir -p data/$PROJECT/output/hsqs/${GENE}
mkdir -p data/$PROJECT/output/mgrms/${GENE}
mkdir -p data/$PROJECT/output/results/${GENE}

I=1
J=1
for ((A=0; A<${SGE_TASK_ID}; A++)); do
    if [[ $J -lt $TRANSCRIPTS ]]; then
        let J=J+1
    else
        let I=I+1
        let J=I+1
    fi
done

L1=$(awk -v I=$I '{if (NR==I) print $1}' data/${PROJECT}/input/${GENE}_ids)
T1=$(awk -v I=$I '{if (NR==I) print $2}' data/${PROJECT}/input/${GENE}_ids)
L2=$(awk -v J=$J '{if (NR==J) print $1}' data/${PROJECT}/input/${GENE}_ids)
T2=$(awk -v J=$J '{if (NR==J) print $2}' data/${PROJECT}/input/${GENE}_ids)
TNAME=$(echo ${T1}_${T2})

# Info
echo pair: $I $J
echo line number: $L1 $L2
echo transcript names: $T1 $T2

LINE=$(sed -n ${L1}p data/${PROJECT}/input/phenotype_ids)
echo $LINE | \
    (
    read ID CHR START END EXTRA
    printf "ID: %s\tCHR: %s\n" $ID $CHR

    # Set window
    W_START=$(( ${START} - ${WINDOW} ))
    W_END=$(( ${END} + ${WINDOW} ))

    if [[ ${W_START} -lt 0 ]]; then
        W_START=0
    fi

    # Create the range file for PLINK to generate the cis and trans regions
    printf "%s %s %s R1" $CHR $W_START $W_END > \
        data/${PROJECT}/output/grm_ranges/${GENE}/$TNAME.txt

    # Create file for GCTA pointing to which PED files belong to the trans region
    printf "data/%s/output/grms/%s/%s_trans_int" $PROJECT $GENE $TNAME > \
        data/${PROJECT}/output/grm_chrs/${GENE}/$TNAME.txt
    for ((i=1;i<=22;i++)); do
        if [[ $i -ne $CHR ]]; then
            printf "\ndata/%s/output/grms/ped_file%s" $PROJECT $i >> \
                data/${PROJECT}/output/grm_chrs/${GENE}/$TNAME.txt
        fi
    done

    # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
    printf "data/%s/output/grms/%s/%s_cis\ndata/%s/output/grms/%s/%s_trans" \
        $PROJECT $GENE $TNAME $PROJECT $GENE $TNAME > \
        data/${PROJECT}/output/mgrms/${GENE}/$TNAME.txt

    # Split the cis and trans regions into separate PED files
    # If no cis SNPs, exit
    if plink --bfile data/${PROJECT}/input/ped_file \
        --extract range data/${PROJECT}/output/grm_ranges/${GENE}/$TNAME.txt \
        --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
        --out data/${PROJECT}/output/fin_peds/${GENE}/${TNAME}_cis; then
        plink --bfile data/${PROJECT}/input/ped_file \
            --exclude range data/${PROJECT}/output/grm_ranges/${GENE}/$TNAME.txt \
            --chr $CHR --remove data/${PROJECT}/input/removed_samples.txt --make-bed \
            --out data/${PROJECT}/output/fin_peds/${GENE}/${TNAME}_trans_int

        # Generate the GRM for the cis and trans region
        gcta64 --make-grm-bin --thread-num $THREADS \
            --bfile data/${PROJECT}/output/fin_peds/${GENE}/${TNAME}_cis \
            --make-grm-alg 0 \
            --out data/${PROJECT}/output/grms/${GENE}/${TNAME}_cis
        gcta64 --make-grm-bin --thread-num $THREADS \
            --bfile data/${PROJECT}/output/fin_peds/${GENE}/${TNAME}_trans_int \
            --make-grm-alg 0 --chr $CHR \
            --out data/${PROJECT}/output/grms/${GENE}/${TNAME}_trans_int
        gcta64 --make-grm-bin --thread-num $THREADS \
            --mgrm data/${PROJECT}/output/grm_chrs/${GENE}/$TNAME.txt \
            --out data/${PROJECT}/output/grms/${GENE}/${TNAME}_trans

    # Run GREML, defaulting to EM if AI fails
    gcta64 --reml-bivar $L1 $L2 --thread-num $THREADS --reml-alg 0 --reml-maxit 10000 \
        --reml-bivar-lrt-rg 0 \
        --remove data/${PROJECT}/input/removed_samples.txt \
        --mgrm data/${PROJECT}/output/mgrms/${GENE}/$TNAME.txt \
        --pheno data/${PROJECT}/input/phenotype \
        --out data/${PROJECT}/output/hsqs/${GENE}/${T1}_${T2} || \
    gcta64 --reml-bivar $L1 $L2 --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 \
        --reml-bivar-lrt-rg 0 \
        --remove data/${PROJECT}/input/removed_samples.txt \
        --mgrm data/${PROJECT}/output/mgrms/${GENE}/$TNAME.txt \
        --pheno data/${PROJECT}/input/phenotype \
        --out data/${PROJECT}/output/hsqs/${GENE}/${T1}_${T2}

        # Add number of cis SNPs to HSQ output
        plink --bfile data/${PROJECT}/input/ped_file \
            --extract range data/${PROJECT}/output/grm_ranges/${GENE}/$TNAME.txt \
            --remove data/${PROJECT}/input/removed_samples.txt --write-snplist \
            --out data/${PROJECT}/output/grm_ranges/${GENE}/$TNAME
        wc -l < data/${PROJECT}/output/grm_ranges/${GENE}/${TNAME}.snplist
        printf "cisSNPs\t%s" $(wc -l < data/${PROJECT}/output/grm_ranges/${GENE}/${TNAME}.snplist) \
            >> data/${PROJECT}/output/hsqs/${GENE}/$TNAME.hsq

        # Clean up files
        rm -rf data/${PROJECT}/output/grms/${GENE}/${TNAME}*
    fi

    # Clean up files
    rm -rf data/${PROJECT}/output/grm_ranges/${GENE}/${TNAME}*
    rm -rf data/${PROJECT}/output/grm_chrs/${GENE}/${TNAME}*
    rm -rf data/${PROJECT}/output/mgrms/${GENE}/${TNAME}*
    rm -rf data/${PROJECT}/output/fin_peds/${GENE}/${TNAME}*
    )

