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
    open(plinkrangefile, "w") do io
        write(io, "$Chr $WindowStart $WindowEnd R1")
    end

    try
        # Generate cis subset of genotype data
        run(`plink \
            --bfile $(snakemake.params["genotypeprefix"]) \
            --extract range $plinkrangefile \
            --remove $(snakemake.input["excludedsamples"]) \
            --out $cisplink \
            --make-bed`)
    catch error
        println("$ID does not have any cis-SNPs")
        touch("$(snakemake.params["outgremlintermediatedir"]).$ID.done")
        return
    end

    # Generate GRMs
    cisgrm = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_cis_grm")
    run(`gcta64 \
        --bfile $cisplink \
        --out $cisgrm \
        --make-grm-bin \
        --make-grm-alg 0`)

    # Run GREML
    gremlmgrmslist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_mgrms_list.txt")
    snplist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_snps")
    open(gremlmgrmslist, "w") do io
        println(io, cisgrm)
    end
    try
        if !success(`gcta64 \
                    --reml \
                    --reml-alg 0 \
                    --reml-maxit 100 \
                    --mpheno $phenotypeind \
                    --mgrm $gremlmgrmslist \
                    --pheno $(snakemake.input["phenotype"]) \
                    --covar $(snakemake.input["quantcov"]) \
                    --qcovar $(snakemake.input["qualcov"]) \
                    --reml-lrt 1 \
                    --out $(snakemake.params["outhsqdir"])$ID`) 
            run(`gcta64 \
                --reml \
                --reml-alg 2 \
                --reml-maxit 10000 \
                --mpheno $phenotypeind \
                --mgrm $gremlmgrmslist \
                --pheno $(snakemake.input["phenotype"]) \
                --covar $(snakemake.input["quantcov"]) \
                --qcovar $(snakemake.input["qualcov"]) \
                --reml-lrt 1 \
                --out $(snakemake.params["outhsqdir"])$ID`)
        end
        run(`plink --bfile $(snakemake.params["genotypeprefix"]) \
            --extract range $plinkrangefile \
            --remove $(snakemake.input["excludedsamples"]) \
            --write-snplist \
            --out $snplist`)
        open("$(snakemake.params["outhsqdir"])$ID.hsq", "a") do io
            write(io, "cisSNPs\t" * readchomp(pipeline(`cat $snplist.snplist`, `wc -l`)))
        end
    catch
        println("GREML for $ID ran into error.")
        return
    end
    println("GREML for $ID completed successfully")
    touch("$(snakemake.params["outgremlintermediatedir"]).$ID.done")
    return
end

function cleanup(phenotypeind::Integer)
    ID, Chr, Start, End, WindowStart, WindowEnd = phenotypeinfodf[phenotypeind,:]
    files = readdir(snakemake.params["outgremlintermediatedir"])
    rmfiles = files[ [startswith(file, ID) for file in files] ]
    for file in rmfiles
        rm(joinpath(snakemake.params["outgremlintermediatedir"], file))
    end
    return
end

phenotypedf, phenotypeinfodf = preparedfs()
for i in 1:size(phenotypeinfodf)[1]
    greml(i)
    cleanup(i)
    flush(stdout) # I think snakemake is buffering the output
end

touch(snakemake.output[1])

