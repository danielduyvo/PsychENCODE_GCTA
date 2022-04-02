#!/usr/bin/env julia

using CSV, DataFrames

# Necessary inputs:
# @param input Dict: keys:
# excludedsamples
# phenotype
# quantcov
# qualcov
# @param params Dict: keys:
# genotypeprefix
# grmprefix
# outgremlintermediatedir
# outhsqdir
# @param phenotypefile AbstractString: Path to phenotype file
# @param phenotypeinfofile AbstractString: Path to phenotype info file
# @param windowsize Integer: Number of base pairs the window extends in
# either direction

function preparedfs(phenotypefile::AbstractString, phenotypeinfofile::AbstractString,
        windowsize::Integer)
    phenotypedf = CSV.read(phenotypefile, DataFrame,
                           header = 0)
    phenotypeinfodf = CSV.read(phenotypeinfofile, DataFrame,
                               header = [:ID, :Chr, :Start, :End])
    transform!(phenotypeinfodf,
               :Start => function(col)
                   newcol = col .- windowsize
                   for i in 1:length(newcol)
                       if newcol[i] < 0
                           newcol[i] = 0
                       end
                   end
                   return(newcol)
               end => :WindowStart,
               :End => function(col)
                   col .+ windowsize
               end => :WindowEnd)
    return(phenotypedf, phenotypeinfodf)
end

# @param input Dict with keys:
# excludedsamples
# phenotype
# quantcov
# qualcov
# @param params Dict with keys:
# genotypeprefix
# grmprefix
# outgremlintermediatedir
# outhsqdir
function greml(
        phenotypeind::Integer,
        input::Dict,
        params::Dict
    )
    ID, Chr, Start, End, WindowStart, WindowEnd = phenotypeinfodf[phenotypeind,:]
    # If GREML has already been done on this phenotype, skip
    if isfile("$(params["outgremlintermediatedir"]).$ID.done")
        println("GREML for $ID has already been performed")
        return
    end
    plinkrangefile = joinpath(params["outgremlintermediatedir"], ID * "plink_range.txt")
    cisplink = joinpath(params["outgremlintermediatedir"], ID * "_cis")
    open(plinkrangefile, "w") do io
        write(io, "$Chr $WindowStart $WindowEnd R1")
    end

    try
        # Generate cis subset of genotype data
        run(`plink \
            --bfile $(params["genotypeprefix"]) \
            --extract range $plinkrangefile \
            --remove $(input["excludedsamples"]) \
            --out $cisplink \
            --make-bed`)
    catch error
        println("$ID does not have any cis-SNPs")
        touch("$(params["outgremlintermediatedir"]).$ID.done")
        return
    end

    # Generate GRMs
    cisgrm = joinpath(params["outgremlintermediatedir"], ID * "_cis_grm")
    run(`gcta64 \
        --bfile $cisplink \
        --out $cisgrm \
        --make-grm-bin \
        --make-grm-alg 0`)

    # Run GREML
    gremlmgrmslist = joinpath(params["outgremlintermediatedir"], ID * "_mgrms_list.txt")
    snplist = joinpath(params["outgremlintermediatedir"], ID * "_snps")
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
                    --pheno $(input["phenotype"]) \
                    --covar $(input["quantcov"]) \
                    --qcovar $(input["qualcov"]) \
                    --reml-lrt 1 \
                    --out $(params["outhsqdir"])$ID`) 
            run(`gcta64 \
                --reml \
                --reml-alg 2 \
                --reml-maxit 10000 \
                --mpheno $phenotypeind \
                --mgrm $gremlmgrmslist \
                --pheno $(input["phenotype"]) \
                --covar $(input["quantcov"]) \
                --qcovar $(input["qualcov"]) \
                --reml-lrt 1 \
                --out $(params["outhsqdir"])$ID`)
        end
        run(`plink --bfile $(params["genotypeprefix"]) \
            --extract range $plinkrangefile \
            --remove $(input["excludedsamples"]) \
            --write-snplist \
            --out $snplist`)
        open("$(params["outhsqdir"])$ID.hsq", "a") do io
            write(io, "cisSNPs\t" * readchomp(pipeline(`cat $snplist.snplist`, `wc -l`)))
        end
    catch
        println("GREML for $ID ran into error.")
        return
    end
    println("GREML for $ID completed successfully")
    touch("$(params["outgremlintermediatedir"]).$ID.done")
    return
end

# @param params includes "outgremlintermediatedir"
function cleanup(phenotypeind::Integer, params::Dict)
    ID, Chr, Start, End, WindowStart, WindowEnd = phenotypeinfodf[phenotypeind,:]
    files = readdir(params["outgremlintermediatedir"])
    rmfiles = files[ [startswith(file, ID) for file in files] ]
    for file in rmfiles
        rm(joinpath(params["outgremlintermediatedir"], file))
    end
    return
end

