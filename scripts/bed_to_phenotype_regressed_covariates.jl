#!/usr/bin/env julia

using CSV, DataFrames, GLM, LinearAlgebra, StatsModels

function generatephenotype_regresscovariates(bedfilename, covfilename)
    beddf = CSV.read(bedfilename, DataFrame)
    chrcolname = names(beddf)[1]
    newchrcolname = replace(chrcolname, r"#" => "")
    rename!(beddf, chrcolname => newchrcolname)
    if !(eltype.(eachcol(beddf))[1] <: Integer)
        subset!(beddf, Symbol(newchrcolname) => chrom -> match(r"\D", chrom) == nothing)
    end

    # FastQTL has 4 columns before the samples, QTLtools has 6
    idcol = 4
    if (eltype.(eachcol(beddf))[5] <: Number)
        firstcol = 5
    elseif (eltype.(eachcol(beddf))[7] <: Number)
        firstcol = 7
    else
        throw("Does not match either FastQTL or QTLTools output")
    end

    cov = CSV.read(covfilename, DataFrame,
                   transpose = true)

    phenotypemat = Matrix(select(beddf, firstcol:size(beddf)[2]))
    phenotypemat = permutedims(phenotypemat, (2, 1)) # matrix is now samples x phenotypes
    cov = cov[indexin(cov.id, names(beddf)[firstcol:end]), :] # just in case
    f = ConstantTerm(0) ~ ConstantTerm(1) + sum([Term(Symbol(n)) for n in names(cov)[2:end]])
    covmat = modelmatrix(f, cov)
    beta = inv(Symmetric(transpose(covmat) * covmat)) * transpose(covmat) * phenotypemat
    phenotypemat .-= covmat * beta
    phenotypedf = DataFrame(phenotypemat, :auto)
    phenotypedf[!, :FamID] = zeros(Integer, size(phenotypedf)[1])
    phenotypedf[!, :SampleID] = names(beddf)[firstcol:size(beddf)[2]]
    select!(phenotypedf, :FamID, :SampleID, :)
    phenotypeinfo = select(beddf, idcol, 1:3)

    return (phenotypedf, phenotypeinfo)
end

phenotypedf1, phenotypeinfo =
generatephenotype_regresscovariates(snakemake.input["pheno1"], snakemake.input["cov1"])
phenotypedf2, _ =
generatephenotype_regresscovariates(snakemake.input["pheno2"], snakemake.input["cov2"])

CSV.write(snakemake.output["pheno"], vcat(phenotypedf1, phenotypedf2),
          delim = '\t',
          writeheader = false)
CSV.write(snakemake.output["pheno_info"], phenotypeinfo,
          delim = '\t',
          writeheader = false)

