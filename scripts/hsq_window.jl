#!/usr/bin/env julia

using CSV,DataFrames

function readhsqs(hsqdir::String)
    files = readdir(hsqdir)
    gremlfiles = files[ [splitext(file)[2] == ".hsq" for file in files] ]
    nresults = length(gremlfiles)
    gremlresults = DataFrame(:ID => Array{Union{Missing, String}}(missing, nresults),
                             :V_Cis => Array{Union{Missing, Float64}}(missing, nresults),
                             :V_Trans => Array{Union{Missing, Float64}}(missing, nresults),
                             :V_Res => Array{Union{Missing, Float64}}(missing, nresults),
                             :V_Phe => Array{Union{Missing, Float64}}(missing, nresults),
                             :Cis_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :Trans_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :SNP_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Cis => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Trans => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Res => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_V_Phe => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_Cis_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_Trans_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :SE_SNP_h2 => Array{Union{Missing, Float64}}(missing, nresults),
                             :P_val => Array{Union{Missing, Float64}}(missing, nresults),
                             :Cis_SNPs => Array{Union{Missing, Int64}}(missing, nresults))
    rown = 1
    for file in gremlfiles
        hsq = CSV.read(joinpath(hsqdir, file), DataFrame, header = 0, delim = '\t', 
                       skipto = 2, types = [String, Float64, Float64],
                       silencewarnings = true)
        hsqrow = (ID = splitext(file)[1]::String,
                  V_Cis = hsq[1,2]::Float64,
                  V_Trans = hsq[2,2]::Float64,
                  V_Res = hsq[3,2]::Float64,
                  V_Phe = hsq[4,2]::Float64,
                  Cis_h2 = hsq[5, 2]::Float64,
                  Trans_h2 = hsq[6, 2]::Float64,
                  SNP_h2 = hsq[7, 2]::Float64,
                  SE_V_Cis = hsq[1,3]::Float64,
                  SE_V_Trans = hsq[2,3]::Float64,
                  SE_V_Res = hsq[3,3]::Float64,
                  SE_V_Phe = hsq[4,3]::Float64,
                  SE_Cis_h2 = hsq[5,3]::Float64,
                  SE_Trans_h2 = hsq[6,3]::Float64,
                  SE_SNP_h2 = hsq[7,3]::Float64,
                  P_val = hsq[12,2]::Float64,
                  Cis_SNPs = hsq[14,2]::Float64)
        gremlresults[rown, :] = hsqrow
        rown = rown + 1
    end
    return gremlresults
end

gremlresultsdf = readhsqs(snakemake.params["outhsqdir"])
CSV.write(snakemake.output[1], gremlresultsdf)
