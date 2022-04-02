#!/usr/bin/env julia

include("scripts/run_cis.jl")

phenotypedf, phenotypeinfodf = preparedfs(snakemake.input["phenotype"],
                                          snakemake.input["phenotypeinfo"],
                                          snakemake.params["windowsize"])
for i in 1:size(phenotypeinfodf)[1]
    greml(i, snakemake.input, snakemake.params)
    cleanup(i, snakemake.params)
    flush(stdout) # I think snakemake is buffering the output
end

touch(snakemake.output[1])

