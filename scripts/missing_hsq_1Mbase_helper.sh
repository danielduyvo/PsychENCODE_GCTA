#!/usr/bin/env bash

# Listing genes missing SNPs in the cis window
> data/${PROJECT}/output/results/missing.txt
for file in $(ls data/${PROJECT}/output/hsqs | grep log | sed 's/\.log/.hsq/g'); do
    if [[ ! -e data/${PROJECT}/output/hsqs/$file ]]; then
        printf '%s\n' $file >> data/${PROJECT}/output/results/missing.txt
    fi
done
