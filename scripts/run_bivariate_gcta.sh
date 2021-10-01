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

while getopts a:t:w:mrdh opt; do
    case $opt in
        h)
            echo 24905 genes
            echo 93293 isoforms
            echo -a for analysis type:
            echo -t for job array string
            echo -r to redo
            echo -m to list missing analyses
            echo -w to set window size
            echo -d to delete generated script file
            echo positional parameters:
            echo POPULATION PHENOTYPE GENOTYPE MODEL
            exit 0
            ;;
        a)
            ANALYSIS=$OPTARG
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
        w)
            window=$OPTARG
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

shift $((OPTIND-1))
POP=$1
PHENO=$2
GENO=$3
MODEL=$4
NAME="${POP}_${PHENO}_${GENO}_${MODEL}"

if $MISSING; then
    echo Option not implemented 1>&2
    exit 1

    echo Generating list of missing HSQs...
    export PROJECT=$NAME
    ./scripts/missing_bivariate_hsq.sh #TODO
fi

if $REDO; then
    ARRAY_STRING=$(printf '1-%s' $(cat data/${NAME}/output/results/missing.txt | wc -l))

    echo Option not implemented 1>&2
    exit 1 #TODO
fi

if [[ -z "$ANALYSIS" || -z "$NAME" || -z "$ARRAY_STRING" || \
    -z "$POP" || -z "$PHENO" || -z "$GENO" || -z "$MODEL" ]]; then
    echo "Missing required option" 1>&2
    exit 1
fi

echo Analysis: $ANALYSIS
echo Project Name: $NAME
echo Job Array: $ARRAY_STRING
echo Redo: $REDO
echo List missing: $MISSING

# Analysis
case $ANALYSIS in
    window)

        ### Read in window size
        if [[ ! -z "$window" ]]; then
            echo Window: $window
        else
            echo "Window option required" 1>&2
            exit 1
        fi

        ## Generate a list of unique groups (genes) of phenotypes (isoforms)
        awk '{print $5}' data/${NAME}/input/phenotype_ids | sort | uniq > \
            data/${NAME}/input/group_ids

        ## Read in number of groups
        unset IFS
        unset group_range
        unset group_indexes
        if [[ -n $(echo $ARRAY_STRING | grep "," | grep "-") ]]; then
            echo "Can only have an array or range, not both" 1>&2
            # exit 1
        elif [[ -n $(echo $ARRAY_STRING | grep "-") ]]; then
            echo "Reading in a range"
            IFS='- '
            read -a group_range <<< $ARRAY_STRING
            unset IFS
            group_indexes=(`seq ${group_range[0]} ${group_range[1]}`)
        else
            echo "Reading in a list of indexes"
            IFS=', '
            read -a group_indexes <<< $ARRAY_STRING
            unset IFS
        fi

        ## Generate qsub script for each group
        for i in "${group_indexes[@]}"; do

            ### Read in group name and create a list of phenotypes
            GROUP_NAME=$(sed -n ${i}p data/${NAME}/input/group_ids)
            sed "s/ARG_POP/${POP}/g; s/ARG_PHENO/${PHENO}/g; s/ARG_GENO/${GENO}/g; s/ARG_MODEL/${MODEL}/g; s/ARG_GENE/${GROUP_NAME}/g" scripts/bivariate_helper.sh > generate_phenotype_group_list.sh
            chmod +x generate_phenotype_group_list.sh
            ./generate_phenotype_group_list.sh

            ### Read in number of phenotypes and convert to array string
            PHENOTYPE_SUM=$(cat data/${NAME}/input/${GROUP_NAME}_ids | wc -l)
            PHENOTYPE_COMBINATIONS=$(( ${PHENOTYPE_SUM} * ( ${PHENOTYPE_SUM} - 1 ) / 2 ))
            if [[ $PHENOTYPE_COMBINATIONS -gt 0 ]]; then
                ARRAY_STRING="1-${PHENOTYPE_COMBINATIONS}"
                echo Group: $GROUP_NAME
                echo Job array: $ARRAY_STRING
                sed "s/ARG_NAME/$NAME/g; s/ARG_ARRAY/$ARRAY_STRING/g; s/ARG_WINDOW/$window/g; s/ARG_GROUP/${GROUP_NAME}/g; s/ARG_REDO/$REDO/g" \
                    scripts/run_bivariate_window.sh > qsub_script.sh
                qsub qsub_script.sh
            fi

        done
        ;;
esac

if $DELETE; then
    rm qsub_script.sh
fi

