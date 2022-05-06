#!/usr/bin/env julia

using CSV, DataFrames, Glob

function readhsqs(hsqdir::String)
    gremlfiles = glob(joinpath(hsqdir, "*.hsq"))
    nresults = length(gremlfiles)
    gremlresults = DataFrame(:ID => Array{String}(undef, nresults),
                             :V_Cis => Array{Float64}(undef, nresults),
                             :V_Res => Array{Float64}(undef, nresults),
                             :V_Phe => Array{Float64}(undef, nresults),
                             :Cis_h2 => Array{Float64}(undef, nresults),
                             :SE_V_Cis => Array{Float64}(undef, nresults),
                             :SE_V_Res => Array{Float64}(undef, nresults),
                             :SE_V_Phe => Array{Float64}(undef, nresults),
                             :SE_Cis_h2 => Array{Float64}(undef, nresults),
                             :P_val => Array{Float64}(undef, nresults),
                             :Cis_SNPs => Array{Int64}(undef, nresults))
    allowmissing!(gremlresults)
    rown = 1
    for file in gremlfiles
        hsq = CSV.read(file, DataFrame; header = 0, delim = '\t', 
                       skipto = 2, types = [String, Float64, Float64],
                       silencewarnings = true)
        hsqrow = (ID = splitext(basename(file))[1]::String,
                  V_Cis = hsq[1,2]::Float64,
                  V_Res = hsq[2,2]::Float64,
                  V_Phe = hsq[3,2]::Float64,
                  Cis_h2 = hsq[4, 2]::Float64,
                  SE_V_Cis = hsq[1,3]::Float64,
                  SE_V_Res = hsq[2,3]::Float64,
                  SE_V_Phe = hsq[3,3]::Float64,
                  SE_Cis_h2 = hsq[4,3]::Float64,
                  P_val = hsq[9,2]::Float64,
                  Cis_SNPs = size(hsq, 1) == 11 ? convert(Int64, hsq[11,2]) : missing)
        gremlresults[rown, :] = hsqrow
        rown = rown + 1
    end
    return gremlresults
end

gremlresultsdf = readhsqs(snakemake.params["outhsqdir"])
CSV.write(snakemake.output[1], gremlresultsdf)
