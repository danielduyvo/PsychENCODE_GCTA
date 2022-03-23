#!/usr/bin/env julia

using CSV,DataFrames

function readhsqs(hsqdir::String)
    files = readdir(hsqdir)
    gremlfiles = files[ [splitext(file)[2] == ".hsq" for file in files] ]
    nresults = length(gremlfiles)
    gremlresults = DataFrame(:ID => Array{Union{Missing, String}}(missing, nresults),
                             :V_Cis => Array{Union{Missing, Float64}}(missing, nresults),
                             :V_Res => Array{Union{Missing, Float64}}(missing, nresults),
                             :V_Phe => Array{Union{Missing, Float64}}(missing, nresults),
                             :Cis_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Cis => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Res => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Phe => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_Cis_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :P_val => Array{Union{Missing, Float64}}(missing, nresults),
                             :Cis_SNPs => Array{Union{Missing, Int64}}(missing, nresults))
    rown = 1
    for file in gremlfiles
        hsq = CSV.read(joinpath(hsqdir, file), DataFrame, header = 0, delim = '\t', 
                       skipto = 2, types = [String, Float64, Float64],
                       silencewarnings = true)
        hsqrow = (ID = splitext(file)[1]::String,
                  V_Cis = hsq[1,2]::Float64,
                  V_Res = hsq[2,2]::Float64,
                  V_Phe = hsq[3,2]::Float64,
                  Cis_h2 = hsq[4, 2]::Float64,
                  SE_V_Cis = hsq[1,3]::Float64,
                  SE_V_Res = hsq[2,3]::Float64,
                  SE_V_Phe = hsq[3,3]::Float64,
                  SE_Cis_h2 = hsq[4,3]::Float64,
                  P_val = hsq[9,2]::Float64,
                  Cis_SNPs = hsq[11,2]::Int64)
        gremlresults[rown, :] = hsqrow
        rown = rown + 1
    end
    return gremlresults
end

gremlresultsdf = readhsqs(snakemake.params["outhsqdir"])
CSV.write(snakemake.output[1], gremlresultsdf)
