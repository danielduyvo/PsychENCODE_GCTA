#!/usr/bin/env bash

POP="EUR"
PHENO="SPC_isoform"
GENO="HRC"
MODEL="bivariate_cis"
PROJECT="${POP}_${PHENO}_${GENO}_${MODEL}"

GENE="ENSG00000128891"

awk -v GENE=$GENE '{if($0 ~ GENE) print NR " " $1}' \
    data/${PROJECT}/input/phenotype_ids > \
    data/${PROJECT}/input/${GENE}_ids

mkdir -p data/${PROJECT}/output/${GENE}/fin_peds
mkdir -p data/${PROJECT}/output/${GENE}/graphs
mkdir -p data/${PROJECT}/output/${GENE}/grm_chrs
mkdir -p data/${PROJECT}/output/${GENE}/grm_ranges
mkdir -p data/${PROJECT}/output/${GENE}/grms
mkdir -p data/${PROJECT}/output/${GENE}/hsqs
mkdir -p data/${PROJECT}/output/${GENE}/mgrms
mkdir -p data/${PROJECT}/output/${GENE}/results

pushd data/$PROJECT/output/${GENE}/grms
ln -s ../../../../../grms/${POP}_${GENO}/* ./
popd
