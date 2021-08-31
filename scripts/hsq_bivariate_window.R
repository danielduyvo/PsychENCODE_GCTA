#!/usr/bin/env Rscript

PROJECT="EUR_SPC_isoform_HRC_bivariate_window"
GENE="ENSG00000128891"

# Combine the HSQ files and return statistics
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/hsqs/", GENE, "/", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                  hsq = read.table(paste('data/', PROJECT, '/output/hsqs/', GENE, '/', file_name, sep=''),
                                   header=TRUE, nrows=17, fill=TRUE)
                  hsq_log = read.table(paste('data/', PROJECT, '/output/hsqs/', GENE, '/', file_name, sep=''),
                                       header=FALSE, skip=18, fill=TRUE)
                  transcripts = substring(file_name, 1, nchar(file_name) - 4)
                  transcripts = unlist(strsplit(transcripts, split="_"))

                  # Transcript ID1, Transcript ID2,
                  # T1 Variance from Cis, T2 Variance from Cis,
                  # Covariance from Cis between T1 and T2,
                  # T1 Variance from Trans, T2 Variance from Trans,
                  # Covariance from Trans between T1 and T2,
                  # T1 Residual Variance, T2 Residual Variance,
                  # Residual Covariance between T1 and T2,
                  # T1 Phenotypic Variance, T2 Phenotypic Variance,
                  # T1 Cis h2, T2 Cis h2,
                  # T1 Trans h2, T2 Trans h2,
                  # Cis Genetic Correlation, Trans Genetic Correlation
                  # SE T1 Variance from Cis, SE T2 Variance from Cis,
                  # SE Covariance from Cis between T1 and T2,
                  # SE T1 Variance from Trans, SE T2 Variance from Trans,
                  # SE Covariance from Trans between T1 and T2,
                  # SE T1 Residual Variance, SE T2 Residual Variance,
                  # SE Residual Covariance between T1 and T2,
                  # SE T1 Phenotypic Variance, SE T2 Phenotypic Variance,
                  # SE T1 Cis h2, SE T2 Cis h2,
                  # SE T1 Trans h2, SE T2 Trans h2,
                  # SE Cis Genetic Correlation, SE Trans Genetic Correlation
                  # P-value
                  line = c(transcripts[1], transcripts[2], 
                           hsq[1,2], hsq[2,2],
                           hsq[3,2],
                           hsq[4,2], hsq[5,2],
                           hsq[6,2],
                           hsq[7,2], hsq[8,2],
                           hsq[9,2],
                           hsq[10,2], hsq[11,2],
                           hsq[12,2], hsq[13,2],
                           hsq[14,2], hsq[15,2],
                           hsq[16,2], hsq[17,2],

                           hsq[1,3], hsq[2,3],
                           hsq[3,3],
                           hsq[4,3], hsq[5,3],
                           hsq[6,3],
                           hsq[7,3], hsq[8,3],
                           hsq[9,3],
                           hsq[10,3], hsq[11,3],
                           hsq[12,3], hsq[13,3],
                           hsq[14,3], hsq[15,3],
                           hsq[16,3], hsq[17,3],

                           hsq_log[5,2]
                  )
                  return(line)
            }))

                  # Transcript ID1, Transcript ID2,
                  # T1 Variance from Cis, T2 Variance from Cis,
                  # Covariance from Cis between T1 and T2,
                  # T1 Variance from Trans, T2 Variance from Trans,
                  # Covariance from Trans between T1 and T2,
                  # T1 Residual Variance, T2 Residual Variance,
                  # Residual Covariance between T1 and T2,
                  # T1 Phenotypic Variance, T2 Phenotypic Variance,
                  # T1 Cis h2, T2 Cis h2,
                  # T1 Trans h2, T2 Trans h2,
                  # Cis Genetic Correlation, Trans Genetic Correlation
                  # SE T1 Variance from Cis, SE T2 Variance from Cis,
                  # SE Covariance from Cis between T1 and T2,
                  # SE T1 Variance from Trans, SE T2 Variance from Trans,
                  # SE Covariance from Trans between T1 and T2,
                  # SE T1 Residual Variance, SE T2 Residual Variance,
                  # SE Residual Covariance between T1 and T2,
                  # SE T1 Phenotypic Variance, SE T2 Phenotypic Variance,
                  # SE T1 Cis h2, SE T2 Cis h2,
                  # SE T1 Trans h2, SE T2 Trans h2,
                  # SE Cis Genetic Correlation, SE Trans Genetic Correlation
                  # P-value

colnames(hsqs) = c("ID_1", "ID_2",
                   "V_Cis_1", "V_Cis_2",
                   "C_Cis",
                   "V_Trans_1", "V_Trans_2",
                   "C_Trans",
                   "V_Res_1", "V_Res_2",
                   "C_Res",
                   "V_Phe_1", "V_Phe_2",
                   "Cis_h2_1", "Cis_h2_2",
                   "Trans_h2_1", "Trans_h2_2",
                   "Cis_rG", "Trans_rG",


                   "SE_V_Cis_1", "SE_V_Cis_2",
                   "SE_C_Cis",
                   "SE_V_Trans_1", "SE_V_Trans_2",
                   "SE_C_Trans",
                   "SE_V_Res_1", "SE_V_Res_2",
                   "SE_C_Res",
                   "SE_V_Phe_1", "SE_V_Phe_2",
                   "SE_Cis_h2_1", "SE_Cis_h2_2",
                   "SE_Trans_h2_1", "SE_Trans_h2_2",
                   "SE_Cis_rG", "SE_Trans_rG",

                   "P_val")
hsqs = as.data.frame(hsqs)
for (i in 3:37) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/", GENE, ".txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)


