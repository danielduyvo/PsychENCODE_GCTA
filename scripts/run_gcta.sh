#!/usr/bin/env bash

# Main script for running GCTA analysis

while getopts a:n:t: opt; do
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

case $ANALYSIS in
    23vc)
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

