#!/usr/bin/env bash

first=false
qtltools=false
THREADS=1

while getopts :fqt:h opt; do
    case $opt in
        h)
            echo Script for generating and populating project directory
            echo Takes in positional arguments to name directories and link files
            echo '$1 population'
            echo '$2 phenotype'
            echo '$3 genotype'
            echo '$4 model'
            echo 'If -f flag used, then input files are generated with additional \
                params'
            echo 'If -q flag used, then bed file is in QTLtools format (default \
                is FastQTL)'
            echo '$5 vcf file'
            echo '$6 bed file'
            echo '$7 number of columns in bed file before samples'
            echo Example: scripts/build_project.sh EUR SPC_gene HRC whole
            exit 0
            ;;
        f)
            first=true
            echo Generating necessary input files for GCTA analysis
            ;;
        q)
            qtltools=true
            ;;
        t)
            THREADS=$OPTARG
            echo Threads: $THREADS
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
pop=$1
pheno=$2
geno=$3
model=$4
project="${pop}_${pheno}_${geno}_${model}"

if [[ -z $model ]]; then
    echo Missing arguments 1>&2
    exit 1
fi

if $first; then
    vcf=$5
    bed=$6
    skip=$7
    if [[ -z $skip ]]; then
        echo Missing arguments 1>&2
        exit 1
    fi
    if $qtltools; then
        echo Reading bed file in QTLtools format
    else
        echo Reading bed file in FastQTL format
    fi

    mkdir -p fin_genotype/${pop}_${geno}/
    mkdir -p grms/${pop}_${geno}/
    mkdir -p fin_RNA_seq/${pop}_${pheno}/

    # Remove nonpolymorphic variants
    bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
        $vcf -o fin_genotype/${pop}_${geno}/subset.vcf

    # Convert VCF to PED
    ## Set family ID all to 0
    plink \
        --vcf fin_genotype/${pop}_${geno}/subset.vcf \
        --recode \
        --out fin_genotype/${pop}_${geno}/ped_file \
        --const-fid --make-bed

    # Generate complete GRM
    > grms/${pop}_${geno}/all_chrs.txt
    for ((i=1;i<=22;i++)); do
        gcta64 --make-grm-bin \
            --bfile fin_genotype/${pop}_${geno}/ped_file \
            --make-grm-alg 0 \
            --chr $i \
            --out grms/${pop}_${geno}/ped_file$i \
            --thread-num ${THREADS} && \
        printf "grms/${pop}_${geno}/ped_file%s\n" $PROJECT $i >> grms/${pop}_${geno}/all_chrs.txt
    done

    gcta64 \
        --mgrm grms/${pop}_${geno}/all_chrs.txt \
        --make-grm-bin \
        --out grms/${pop}_${geno}/complete

    # Generate Big K GRM
    gcta64 --grm grms/${pop}_${geno}/complete \
        --thread-num $THREADS \
        --make-bK 0.025 \
        --out grms/${pop}_${geno}/bigk

    # Subset GRM to filter out cryptic relatedness (GR > 0.05)
    gcta64 \
        --grm grms/${pop}_${geno}/complete \
        --grm-cutoff 0.05 \
        --make-grm-bin \
        --out grms/${pop}_${geno}/unrelated
    diff \
        --new-line-format=%L \
        --unchanged-line-format="" \
        grms/${pop}_${geno}/unrelated.grm.id \
        grms/${pop}_${geno}/complete.grm.id |
    awk '{print $1 " " $2}' > \
    fin_genotype/${pop}_${geno}/removed_samples.txt

    # Generate phenotype file
    if $qtltools; then
        Rscript scripts/generate_phenotype_QTLtools.R $bed $pop $pheno
    else
        Rscript scripts/generate_phenotype_FastQTL.R $bed $pop $pheno
    fi

    echo Done generating input files
fi

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
