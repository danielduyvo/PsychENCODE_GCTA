#!/usr/bin/env julia
using JSON

snakemake = JSON.parse(ARGS[2])
include("scripts/run_cis.jl")

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

