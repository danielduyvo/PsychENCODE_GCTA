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
MISSING=false

while getopts a:n:t:w:mrdh opt; do
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
        m)
            MISSING=true
            ;;
        h)
            echo 24905 genes
            echo 93293 isoforms
            echo -a for analysis type:
            echo 23vc, 5vc, bksk, chr, cis, whole, window, window_exc
            echo "-n for project name"
            echo -t for job array string
            echo -r to redo
            echo -m to list missing analyses
            echo -w to set window size
            echo -d to delete generated script file
            exit 0
            ;;
        w)
            WINDOW=$OPTARG
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

if $REDO; then
    ARRAY_STRING=$(printf '1-%s' $(echo data/${NAME}/output/results/missing.txt | wc -l))
fi

if [[ -z "$ANALYSIS" || -z "$NAME" || -z "$ARRAY_STRING" ]]; then
    echo "Missing required option" 1>&2
    exit 1
fi

echo Analysis: $ANALYSIS
echo Project Name: $NAME
echo Job Array: $ARRAY_STRING
echo Redo: $REDO
echo List missing: $MISSING

if $MISSING; then
    export PROJECT
    ./scripts/missing_hsq.sh
fi


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
        sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_REDO/$REDO/g" \
            scripts/run_chr.sh > qsub_script.sh
        qsub qsub_script.sh
        ;;
    cis)
        if [[ ! -z "$WINDOW" ]]; then
            echo WINDOW: $WINDOW
        else
            echo "Window option required" 1>&2
            exit 1
        fi
        sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_WINDOW/$WINDOW/g; s/ARG_REDO/$REDO/g" \
            scripts/run_cis.sh > qsub_script.sh
        qsub qsub_script.sh
        ;;
    whole)
        sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_REDO/$REDO/g" \
            scripts/run_whole.sh > qsub_script.sh
        qsub qsub_script.sh
        ;;
    window)
        if [[ ! -z "$WINDOW" ]]; then
            echo WINDOW: $WINDOW
        else
            echo "Window option required" 1>&2
            exit 1
        fi
        sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_WINDOW/$WINDOW/g; s/ARG_REDO/$REDO/g" \
            scripts/run_window.sh > qsub_script.sh
        qsub qsub_script.sh
        ;;
    window_exc)
        if [[ ! -z "$WINDOW" ]]; then
            echo WINDOW: $WINDOW
        else
            echo "Window option required" 1>&2
            exit 1
        fi
        sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_WINDOW/$WINDOW/g; s/ARG_REDO/$REDO/g" \
            scripts/run_window_cis_chr_excluded_from_trans.sh > qsub_script.sh
        qsub qsub_script.sh
        ;;
esac

if $DELETE; then
    rm qsub_script.sh
fi
