#!/usr/bin/env Rscript

PROJECT="EUR_SPC_isoform_HRC_bivariate_whole"
GENE="ENSG00000078328"

# Combine the HSQ files and return statistics
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/", GENE, "/hsqs", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                  hsq = read.table(paste('data/', PROJECT, '/output/', GENE, '/hsqs/', file_name, sep=''),
                                   header=TRUE, nrows=11, fill=TRUE)
                  hsq_log = read.table(paste('data/', PROJECT, '/output/', GENE, '/hsqs/', file_name, sep=''),
                                       header=FALSE, skip=12, fill=TRUE)
                  transcripts = substring(file_name, 1, nchar(file_name) - 4)
                  transcripts = unlist(strsplit(transcripts, split="_"))

                  # Transcript ID1, Transcript ID2,
                  # T1 Variance from SNP, T2 Variance from SNP,
                  # Covariance from SNP between T1 and T2,
                  # T1 Residual Variance, T2 Residual Variance,
                  # Residual Covariance between T1 and T2,
                  # T1 Phenotypic Variance, T2 Phenotypic Variance,
                  # T1 SNP h2, T2 SNP h2, Genetic Correlation,
                  # SE T1 Variance from SNP, SE T2 Variance from SNP,
                  # SE Covariance from SNP between T1 and T2,
                  # SE T1 Residual Variance, SE T2 Residual Variance,
                  # SE Residual Covariance between T1 and T2,
                  # SE T1 Phenotypic Variance, SE T2 Phenotypic Variance,
                  # SE T1 SNP h2, SE T2 SNP h2, SE Genetic Correlation,
                  # P-value
                  line = c(transcripts[1], transcripts[2], 
                           hsq[1,2], hsq[2,2],
                           hsq[3,2],
                           hsq[4,2], hsq[5,2],
                           hsq[6,2],
                           hsq[7,2], hsq[8,2],
                           hsq[9,2], hsq[10,2], hsq[11,2],

                           hsq[1,3], hsq[2,3],
                           hsq[3,3],
                           hsq[4,3], hsq[5,3],
                           hsq[6,3],
                           hsq[7,3], hsq[8,3],
                           hsq[9,3], hsq[10,3], hsq[11,3],

                           hsq_log[5,2]
                  )
                  return(line)
            }))
                  # Transcript ID1, Transcript ID2,
                  # T1 Variance from SNP, T2 Variance from SNP,
                  # Covariance from SNP between T1 and T2,
                  # T1 Residual Variance, T2 Residual Variance,
                  # Residual Covariance between T1 and T2,
                  # T1 Phenotypic Variance, T2 Phenotypic Variance
                  # T1 SNP h2, T2 SNP h2, Genetic Correlation
                  # SE T1 Variance from SNP, SE T2 Variance from SNP,
                  # SE Covariance from SNP between T1 and T2,
                  # SE T1 Residual Variance, SE T2 Residual Variance,
                  # SE Residual Covariance between T1 and T2,
                  # SE T1 Phenotypic Variance, SE T2 Phenotypic Variance
                  # SE T1 SNP h2, SE T2 SNP h2, SE Genetic Correlation
                  # P-value
colnames(hsqs) = c("ID_1", "ID_2",
                   "V_SNP_1", "V_SNP_2",
                   "C_SNP",
                   "V_Res_1", "V_Res_2",
                   "C_Res",
                   "V_Phe_1", "V_Phe_2",
                   "SNP_h2_1", "SNP_h2_2", "rG",

                   "SE_V_SNP_1", "SE_V_SNP_2",
                   "SE_C_SNP",
                   "SE_V_Res_1", "SE_V_Res_2",
                   "SE_C_Res",
                   "SE_V_Phe_1", "SE_V_Phe_2",
                   "SE_SNP_h2_1", "SE_SNP_h2_2", "SE_rG",

                   "P_val")
hsqs = as.data.frame(hsqs)
for (i in 3:25) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/", GENE, "/results/all_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

