#!/usr/bin/env bash

# Script for compiling GCTA results

while getopts a:n:t:w:mrdh opt; do
    case $opt in
        h)
            echo -a for analysis type:
            echo 23vc, 5vc, bksk, chr, cis, whole, window, window_exc
            echo "-n for project name"
            echo "-w for whole genome all SNPs project name (if applicable)"
            echo -d to delete generated script file
            exit 0
            ;;
        a)
            analysis=$OPTARG
            ;;
        n)
            name=$OPTARG
            ;;
        w)
            whole=$OPTARG
            ;;
        d)
            delete=true
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

if [[ -z "$analysis" || -z "$name" ]]; then
    echo "Missing required option" 1>&2
    exit 1
fi

echo Analysis: $analysis
echo Project Name: $name
echo Delete: $delete

case $analysis in
    23vc)
        sed "s/ARG_NAME/$name/g" \
            scripts/hsq_23vc.R > hsq_compile.R
        Rscript hsq_compile.R
        ;;
    5vc)
        ;;
    bksk)
        ;;
    chr)
        sed "s/ARG_NAME/$name/g" \
            scripts/hsq_chr.R > hsq_compile.R
        Rscript hsq_compile.R
        ;;
    cis)
        sed "s/ARG_NAME/$name/g" \
            scripts/hsq_cis.R > hsq_compile.R
        Rscript hsq_compile.R
        ;;
    whole)
        sed "s/ARG_NAME/$name/g" \
            scripts/hsq_whole.R > hsq_compile.R
        Rscript hsq_compile.R
        ;;
    window)
        if [[ ! -z "$whole" ]]; then
            echo Window: $whole
        else
            echo "Window option required" 1>&2
            exit 1
        fi
        sed "s/ARG_NAME/$name/g; s/ARG_WHOLE/$whole/g" \
            scripts/hsq_window.R > hsq_compile.R
        Rscript hsq_compile.R
        ;;
    window_exc)
        if [[ ! -z "$whole" ]]; then
            echo Window: $whole
        else
            echo "Window option required" 1>&2
            exit 1
        fi
        sed "s/ARG_NAME/$name/g; s/ARG_WHOLE/$whole/g" \
            scripts/hsq_window.R > hsq_compile.R
        Rscript hsq_compile.R
        ;;
esac

if $DELETE; then
    rm hsq_compile.R
fi

