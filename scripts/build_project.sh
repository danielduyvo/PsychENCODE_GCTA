#!/usr/bin/env bash

while getopts :h opt; do
    case $opt in
        h)
            echo Script for generating and populating project directory
            echo Takes in positional arguments to name directories and link files
            echo $1 population
            echo $2 phenotype
            echo $3 genotype
            echo $4 model
            echo Example: scripts/build_project.sh EUR SPC_gene HRC whole
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            exit 1
            ;;
    esac
done

pop=$1
pheno=$2
geno=$3
model=$4

if [[ -z $model ]]; then
    echo Missing arguments 1>&2
    exit 1
fi
project="${pop}_${pheno}_${geno}_${model}"
threads=2
mkdir -p data/$project/output/fin_peds
mkdir -p data/$project/output/graphs
mkdir -p data/$project/output/grm_chrs
mkdir -p data/$project/output/grm_ranges
mkdir -p data/$project/output/grms
mkdir -p data/$project/output/hsqs
mkdir -p data/$project/output/mgrms
mkdir -p data/$project/output/results
mkdir -p data/$project/input/

pushd data/$project/input/
ln -s ../../../fin_genotype/${pop}_${geno}/* ./
ln -s ../../../fin_RNA_seq/${pop}_${pheno}/* ./
popd

pushd data/$project/output/grms
ln -s ../../../../grms/${pop}_${geno}/* ./
popd

echo $project
