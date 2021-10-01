#!/usr/bin/env bash

POP="ARG_POP"
PHENO="ARG_PHENO"
GENO="ARG_GENO"
MODEL="ARG_MODEL"
PROJECT="${POP}_${PHENO}_${GENO}_${MODEL}"

GENE="ARG_GENE"

awk -v GENE=$GENE '{if($0 ~ GENE) print NR " " $1}' \
    data/${PROJECT}/input/phenotype_ids > \
    data/${PROJECT}/input/${GENE}_ids

mkdir -p data/${PROJECT}/output/fin_peds/${GENE}/
mkdir -p data/${PROJECT}/output/graphs/${GENE}/
mkdir -p data/${PROJECT}/output/grm_chrs/${GENE}/
mkdir -p data/${PROJECT}/output/grm_ranges/${GENE}/
mkdir -p data/${PROJECT}/output/grms/${GENE}/
mkdir -p data/${PROJECT}/output/hsqs/${GENE}/
mkdir -p data/${PROJECT}/output/mgrms/${GENE}/
mkdir -p data/${PROJECT}/output/results/${GENE}/

pushd data/$PROJECT/output/grms/${GENE}
ln -s ../../../../grms/${POP}_${GENO}/* ./
popd
