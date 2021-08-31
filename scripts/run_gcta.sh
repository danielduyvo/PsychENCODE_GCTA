#!/usr/bin/env bash

# Main script for running GCTA analysis

# Common -t arguments
# 1-24905
# 1-75000
# 75001-93293

# 24905 genes
# 93293 isoforms

REDO=false

while getopts a:n:t:r opt; do
    case $opt in
        a)
            ANALYSIS=$OPTARG
            ;;
        n)
            NAME=$OPTARG
            ;;
        t)
            ARRAY_STRING=$OPTARG
            ;;
        r)
            REDO=true
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

echo Analysis: $ANALYSIS
echo Project Name: $NAME
echo Job Array: $ARRAY_STRING
echo Redo: $REDO

case $ANALYSIS in
    23vc)
        qsub <(sed "s/ARG_NAME/$NAME; s/ARG_ARRAY/$ARRAY_STRING; s/ARG_REDO/$REDO" \
            scripts/run_23vc.sh)
        ;;
    5vc)
        ;;
    bksk)
        ;;
    chr)
        ;;
    cis)
        ;;
    whole)
        ;;
    window)
        ;;
esac

