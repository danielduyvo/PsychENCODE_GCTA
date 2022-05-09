#!/usr/bin/env julia

using CSV, DataFrames, LinearAlgebra

function readgrmbin(prefix; alln = true, numtype = Float32)
    function sumi(i)::Integer
        return (i + 1) * i / 2
    end
    binfilename = prefix * ".grm.bin"
    nfilename = prefix * ".grm.N.bin"
    idfilename = prefix * ".grm.id"
    id = CSV.read(idfilename, DataFrame;
                  header = 0)
    n = size(id, 1)
    grm = Vector{numtype}(undef, sumi(n))
    read!(binfilename, grm)
    N = Vector{numtype}(undef, sumi(n))
    if alln
        read!(nfilename, N)
    else
        N = Vector{numtype}(undef, 1)
        read!(nfilename, N)
    end
    return (grm, N, id)
end

function writegrmbin(prefix, grm, N, id)
    function sumi(i)::Integer
        return (i + 1) * i / 2
    end
    binfilename = prefix * ".grm.bin"
    nfilename = prefix * ".grm.N.bin"
    idfilename = prefix * ".grm.id"
    id = CSV.write(idfilename, id;
                   delim = "\t", writeheader = false)
    write(binfilename, grm)
    write(nfilename, N)
    return nothing
end

function vectolowertri(vec)
    lvec = length(vec)
    n = convert(Int, 1 / 2 * (sqrt(8 * lvec + 1) - 1))
    @assert (n + 1) * n / 2 == lvec
    i = 1
    j = 1
    mat = Matrix{eltype(vec)}(undef, n, n)
    for el in vec
        mat[i, j] = el
        if j + 1 > i
            i += 1
            j = 1
        else
            j += 1
        end
    end
    return LowerTriangular(mat)
end

function lowertritovec(mat)
    @assert size(mat, 1) == size(mat, 2)
    n = size(mat, 1)
    vec = Vector{eltype(mat)}(undef, convert(Int, (n + 1) * n / 2))
    curr = 1
    for i in 1:n
        for j in 1:i
            vec[curr] = mat[i, j]
            curr += 1
        end
    end
    return vec
end

function generatetri1_tri1tri2(prefix, idfilename)
    grm, N, id = readgrmbin(prefix)
    grmmat = vectolowertri(grm)

    # Generate sample trimester vector
    tri1ids = filter(id -> id != "", split(read(idfilename, String), "\n"))
    sampletrimester = fill(2, size(id, 1))
    for i in 1:length(sampletrimester)
        if id[i, :Column2] in tri1ids
            sampletrimester[i] = 1
        end
    end

    # Edit grm matrices
    tri1_grmmat = deepcopy(grmmat)
    for i in 1:size(tri1_grmmat, 1)
        for j in 1:i
            if sampletrimester[i] == 2 || sampletrimester[j] == 2
                tri1_grmmat[i, j] = 0
            end
        end
    end
    tri1tri2_grmmat = deepcopy(grmmat)
    for i in 1:size(tri1tri2_grmmat, 1)
        for j in 1:i
            if sampletrimester[i] != sampletrimester[j]
                tri1tri2_grmmat[i, j] = 0
            end
        end
    end
    # Convert back to grm vector and print out
    tri1_grmvec = lowertritovec(tri1_grmmat)
    tri1tri2_grmvec = lowertritovec(tri1tri2_grmmat)
    writegrmbin(prefix * "_tri1", tri1_grmvec, N, id)
    writegrmbin(prefix * "_tri1tri2", tri1tri2_grmvec, N, id)
    return nothing
end

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
    generatetri1_tri1tri2(cisgrm, snakemake.input["trisamples"])
    tri1_grm = joinpath(snakemake.params["outgremlintermediatedir"],
                        ID * "_cis_grm" * "_tri1")
    tri1tri2_grm = joinpath(snakemake.params["outgremlintermediatedir"],
                            ID * "_cis_grm" * "_tri1tri2")

    # Run GREML
    gremlmgrmslist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_mgrms_list.txt")
    snplist = joinpath(snakemake.params["outgremlintermediatedir"], ID * "_snps")
    open(gremlmgrmslist, "w") do io
        println(io, tri1_grm)
        println(io, tri1tri2_grm)
    end
    try
        if !success(`gcta64 \
                    --reml \
                    --reml-alg 0 \
                    --reml-maxit 100 \
                    --mpheno $phenotypeind \
                    --mgrm $gremlmgrmslist \
                    --pheno $(snakemake.input["phenotype"]) \
                    --reml-lrt 1 2 \
                    --out $(snakemake.params["outhsqdir"])$ID`) 
            run(`gcta64 \
                --reml \
                --reml-alg 2 \
                --reml-maxit 10000 \
                --mpheno $phenotypeind \
                --mgrm $gremlmgrmslist \
                --pheno $(snakemake.input["phenotype"]) \
                --reml-lrt 1 2 \
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

