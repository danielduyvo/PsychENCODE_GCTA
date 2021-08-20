#!/usr/bin/env bash

# Listing genes missing SNPs in the cis window
> data/${PROJECT}/output/results/missing.txt
for file in $(awk '{print $1}' data/${PROJECT}/input/phenotype_ids); do
    if [[ ! -e data/${PROJECT}/output/hsqs/${file}.hsq ]]; then
        printf '%s.hsq\n' $file >> data/${PROJECT}/output/results/missing.txt
    fi
done
