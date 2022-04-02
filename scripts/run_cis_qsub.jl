#!/usr/bin/env julia
using JSON

snakemake = JSON.parsefile(ARGS[2])
struct Snakemake
    input::Dict
    params::Dict
    output::AbstractVector
end
snakemake = Snakemake(snakemake["input"], snakemake["params"], snakemake["output"])
include("run_cis.jl")

phenotypedf, phenotypeinfodf = preparedfs(snakemake.input["phenotype"],
                                          snakemake.input["phenotypeinfo"
                                                         ][parse(Int64, ARGS[1])],
                                          snakemake.params["windowsize"])
for i in 1:size(phenotypeinfodf)[1]
    greml(i, snakemake.input, snakemake.params)
    cleanup(i, snakemake.params)
    flush(stdout) # I think snakemake is buffering the output
end

touch(snakemake.output[parse(Int64, ARGS[1])])

