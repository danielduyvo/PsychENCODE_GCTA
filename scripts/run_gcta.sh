#!/usr/bin/env bash

# Main script for running GCTA analysis

# Common -t arguments
# 1-24905
# 1-75000
# 75001-93293

# 24905 genes
# 93293 isoforms

REDO=false
DELETE=false

while getopts a:n:t:rdh opt; do
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
        d)
            DELETE=true
            ;;
        h)
            echo 24905 genes
            echo 93293 isoforms
            echo -a for analysis type:
            echo 23vc, 5vc, bksk, chr, cis, whole, window
            echo "-n for project name"
            echo -t for job array string
            echo -r to redo
            echo -d to delete generated script file
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

echo Analysis: $ANALYSIS
echo Project Name: $NAME
echo Job Array: $ARRAY_STRING
echo Redo: $REDO

case $ANALYSIS in
    23vc)
        sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_REDO/$REDO/g" \
            scripts/run_23vc.sh > qsub_script.sh
        qsub qsub_script.sh
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

if $DELETE; then
    rm qsub_script.sh
fi
