#!/usr/env/bin julia

using CSV,DataFrames

processcov = function(covfilename, qualfilename, quantfilename)
    cov = CSV.read(covfilename, DataFrame,
                   transpose = true)
    cov[:,:famid] .= 0
    qualcov = cov[:, unique(vcat("famid", "id",
                                  names(cov, AbstractString)
                                 ))
                  ]
    quantcov = cov[:, unique(vcat("famid", "id",
                                   names(cov, Number)
                                  ))
                   ]

    CSV.write(qualfilename,
              qualcov;
              delim='\t',
              header=false
             )

    CSV.write(quantfilename,
              quantcov;
              delim='\t',
              header=false
             )
    return
end

processcov(snakemake.input[1], snakemake.output[1], snakemake.output[2])

# processcov(
#             "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/chuanjiao/tri1_25HCP.txt",
#             "fin_RNA_seq/Trimester_1_EUR_50_HPC_gene")
# 
# processcov(
#             "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/chuanjiao/tri2_15HCP.txt",
#             "fin_RNA_seq/Trimester_2_EUR_50_HPC_gene")
