#!/bin/bash

# qsub options
#$ -w e
#$ -N EUR_SPC_isoform_HRC_bivariate
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o EUR_SPC_isoform_HRC_bivariate.log
#$ -e EUR_SPC_isoform_HRC_bivariate.err
#$ -m a
#$ -M danieldu
#$ -t 1-28

# Requires a complete GRM to be generated at data/${PROJECT}/output/grms/complete
# Can be generated with the following line using the file data/${PROJECT}/output/grm_chrs/all_chrs.txt
# gcta64 --make-grm-bin --thread-num $THREADS --mgrm data/${PROJECT}/output/grm_chrs/all_chrs.txt \
#     --out data/${PROJECT}/output/grms/complete

PROJECT="EUR_SPC_isoform_HRC_bivariate_whole"
GENE="ENSG00000128891"
TRANSCRIPTS=$(wc -l < data/${PROJECT}/input/${GENE}_ids)
THREADS=2

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

echo pair: $I $J \n
echo line number: $L1 $L2 \n
echo transcript names: $T1 $T2 \n
# Run GREML, defaulting to EM if AI fails
gcta64 --reml-bivar $L1 $L2 --thread-num $THREADS --reml-alg 0 --reml-maxit 10000 \
    --reml-bivar-lrt-rg 0 \
    --remove data/${PROJECT}/input/removed_samples.txt \
    --grm data/${PROJECT}/output/${GENE}/grms/complete\
    --pheno data/${PROJECT}/input/phenotype \
    --out data/${PROJECT}/output/${GENE}/hsqs/${T1}_${T2} || \
gcta64 --reml-bivar $L1 $L2 --thread-num $THREADS --reml-alg 2 --reml-maxit 10000 \
    --reml-bivar-lrt-rg 0 \
    --remove data/${PROJECT}/input/removed_samples.txt \
    --grm data/${PROJECT}/output/${GENE}/grms/complete\
    --pheno data/${PROJECT}/input/phenotype \
    --out data/${PROJECT}/output/${GENE}/hsqs/${T1}_${T2}
