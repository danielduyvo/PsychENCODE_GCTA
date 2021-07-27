#!/bin/bash
set -x
cat test_phenotype_ids | while read line; do
echo $line | \
    (
        read ID CHR START END W_START W_END
        printf "ID: %s\tCHR: %s\n" $ID $CHR

        # Check if the chromosome is a sex chromosome
        if [[ $CHR -eq 0 ]]; then
            continue
        fi

        # Create the range file for PLINK to generate the cis and trans regions
        printf "%s %s %s R1" $CHR $W_START $W_END > grm_ranges/$ID.txt

        # Create file for GCTA pointing to which PED files belong to the trans region
        printf "grms/%s_trans_int" $ID > grm_chrs/$ID.txt
        for ((i=1;i<=22;i++)); do
            if [[ $i -ne $CHR ]]; then
                printf "\ngrms/ped_file%s" $i >> grm_chrs/$ID.txt
            fi
        done

        # Create file for GCTA pointing to which GRM files to compare (cis vs. trans)
        printf "grms/%s_cis\ngrms/%s_trans" $ID $ID > mgrms/$ID.txt

        # Split the cis and trans regions into separate PED files
        plink --bfile ped_file --extract range grm_ranges/$ID.txt --remove removed_samples.txt --make-bed --out fin_peds/${ID}_cis
        plink --bfile ped_file --exclude range grm_ranges/$ID.txt --chr $CHR --remove removed_samples.txt --make-bed --out fin_peds/${ID}_trans_int

        # Generate the GRMs for the cis and trans regions
        gcta64 --make-grm-bin --bfile fin_peds/${ID}_cis --make-grm-alg 0 --out grms/${ID}_cis
        gcta64 --make-grm-bin --bfile fin_peds/${ID}_trans_int --make-grm-alg 0 --chr $CHR --out grms/${ID}_trans_int
        gcta64 --make-grm-bin --mgrm grm_chrs/$ID.txt --out grms/${ID}_trans

        # Run GREML, defaulting to EM if AI fails
        gcta64 --reml --reml-alg 0 --reml-maxit 100 --mpheno 1 --mgrm mgrms/$ID.txt --pheno phenotype --out hsqs/$ID || \
            gcta64 --reml --reml-alg 2 --reml-maxit 10000 --mpheno 1 --mgrm mgrms/$ID.txt --pheno phenotype --out hsqs/$ID

        # Clean up files
        # rm grm_ranges/$ID.txt
        # rm grm_chrs/$ID.txt
        # rm mgrms/$ID.txt
        # rm fin_peds/${ID}_*
        # rm grms/${ID}_*
    )
done
