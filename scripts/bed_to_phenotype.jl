#!/usr/bin/env julia

using CSV,DataFrames

generatephenotypefiles = function(bedfilename, phenotypefilename, infofilename)
    beddf = CSV.read(bedfilename, DataFrame)
    chrcolname = names(beddf)[1]
    newchrcolname = replace(chrcolname, r"#" => "")
    rename!(beddf, chrcolname => newchrcolname)
    if !(Int <: eltype.(eachcol(beddf))[1])
        subset!(beddf, Symbol(newchrcolname) => chrom -> match(r"\D", chrom) == nothing)
    end

    # FastQTL has 4 columns before the samples, QTLtools has 6
    idcol = 4
    if (Number <: eltype.(eachcol(beddf))[5])
        firstcol = 5
    else
        firstcol = 7
    end

    phenotypemat = Matrix(select(beddf, firstcol:size(beddf)[2]))
    phenotypedf = DataFrame(permutedims(phenotypemat, (2, 1)), :auto)
    phenotypedf[!, :FamID] = zeros(Integer, size(phenotypedf)[1])
    phenotypedf[!, :SampleID] = names(beddf)[firstcol:size(beddf)[2]]
    select!(phenotypedf, :FamID, :SampleID, :)
    phenotypeinfo = select(beddf, idcol, 1:3)

    CSV.write(phenotypefilename, phenotypedf,
              delim = '\t',
              writeheader = false)
    CSV.write(infofilename, phenotypeinfo,
              delim = '\t',
              writeheader = false)
end

generatephenotypefiles(snakemake.input[1], snakemake.output[1], snakemake.output[2])
