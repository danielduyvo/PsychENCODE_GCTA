#!/usr/bin/env julia

using CSV,DataFrames

function preparedfs()
    phenotypedf = CSV.read(snakemake.input["phenotype"], DataFrame,
                           header = 0)
    phenotypeinfodf = CSV.read(snakemake.input["phenotypeinfo"], DataFrame,
                               header = [:ID, :Chr, :Start, :End])
    transform!(phenotypeinfodf,
               :Start => function(col)
                   newcol = col .- snakemake.params["windowsize"]
                   for i in 1:length(newcol)
                       if newcol[i] < 0
                           newcol[i] = 0
                       end
                   end
                   return(newcol)
               end => :WindowStart,
               :End => function(col)
                   col .+ snakemake.params["windowsize"]
               end => :WindowEnd)
    return(phenotypedf, phenotypeinfodf)
end

function greml(phenotypeind::Integer)
    ID, Chr, Start, End, WindowStart, WindowEnd = phenotypeinfodf[phenotypeind,:]
    # If GREML has already been done on this phenotype, skip
    if isfile("$(snakemake.params["outgremlintermediatedir"]).$ID.done")
        println("GREML for $ID has already been performed")
        return
    end
    plinkrangefile = joinpath(snakemake.params["outgremlintermediatedir"], ID * "plink_range.txt")
    cisplink = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_cis")
    transintplink = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_trans_int")
    open(plinkrangefile, "w") do io
        write(io, "$Chr $WindowStart $WindowEnd R1")
    end

    try
        # Generate trans and cis subsets of genotype data
        run(`plink \
            --bfile $(snakemake.params["genotypeprefix"]) \
            --extract range $plinkrangefile \
            --remove $(snakemake.input["excludedsamples"]) \
            --out $cisplink \
            --make-bed`)
    catch error
        print("$ID does not have any cis-SNPs")
        touch("$(snakemake.params["outgremlintermediatedir"]).$ID.done")
        return
    end

    run(`plink \
        --bfile $(snakemake.params["genotypeprefix"]) \
        --exclude range $plinkrangefile \
        --chr $Chr \
        --remove $(snakemake.input["excludedsamples"]) \
        --out $transintplink \
        --make-bed`)

    # Generate GRMs
    transmgrmslist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_trans_mgrms_list.txt")
    cisgrm = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_cis_grm")
    transintgrm = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_trans_int_grm")
    transgrm = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_trans_grm")
    open(transmgrmslist, "w") do io
        println(io, transintgrm)
        for i in 1:22
            println(io, "$(snakemake.params["grmprefix"])$i")
        end
    end
    run(`gcta64 \
        --bfile $cisplink \
        --out $cisgrm \
        --make-grm-bin \
        --make-grm-alg 0`)
    run(`gcta64 \
        --bfile $transintplink \
        --out $transintgrm \
        --chr $Chr \
        --make-grm-bin \
        --make-grm-alg 0`)
    run(`gcta64 \
        --mgrm $transmgrmslist \
        --out $transgrm \
        --make-grm-bin`)

    # Run GREML
    gremlmgrmslist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_mgrms_list.txt")
    snplist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "snps")
    open(gremlmgrmslist, "w") do io
        println(io, cisgrm)
        println(io, transgrm)
    end
    if success(`gcta64 \
               --reml \
               --reml-alg 0 \
               --reml-maxit 10000 \
               --mpheno $phenotypeind \
               --mgrm $gremlmgrmslist \
               --pheno $(snakemake.input["phenotype"]) \
               --covar $(snakemake.input["qualcov"]) \
               --qcovar $(snakemake.input["quantcov"]) \
               --reml-lrt 1 2 \
               --out $(snakemake.params["outhsqdir"])$ID "||" \
               gcta64 \
               --reml \
               --reml-alg 2 \
               --reml-maxit 10000 \
               --mpheno $phenotypeind \
               --mgrm $gremlmgrmslist \
               --pheno $(snakemake.input["phenotype"]) \
               --covar $(snakemake.input["qualcov"]) \
               --qcovar $(snakemake.input["quantcov"]) \
               --reml-lrt 1 2 \
               --out $(snakemake.params["outhsqdir"])$ID`)
        run(`plink --bfile $(snakemake.params["genotypeprefix"]) \
                        --extract range $plinkrangefile \
                        --remove $(snakemake.input["excludedsamples"]) \
                        --write-snplist \
                        --out $snplist`)
        open(snplist, "a") do io
            write(io, "cisSNPs\t" * readchomp(`wc -l '<' $snplist`))
        end
    else
        println("GREML for $ID ran into error.")
    end
    println("GREML for $ID completed successfully")
    touch("$(snakemake.params["outgremlintermediatedir"]).$ID.done")
    return
end

phenotypedf, phenotypeinfodf = preparedfs()
for i in 1:size(phenotypeinfodf)[1]
    greml(i)
end
touch("$(snakemake.params["outgremlintermediatedir"]).$(snakemake.wildcards["chunk"])_chunk.done")
