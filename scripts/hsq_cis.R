#!/usr/bin/env Rscript

PROJECT="EUR_SPC_HRC_cis_1Mbase"

# Combine the HSQ files and return statistics
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/hsqs", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                  hsq = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=TRUE, nrows=4, fill=TRUE)
                  hsq_log = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=FALSE, skip=5, fill=TRUE)
                  # Gene ID, 
                  # Variance from Cis, Residual Variance, Phenotypic Variance, Cish2,
                  # SE(Variance from Cis), SE(Residual Variance), SE(Phenotypic Variance), SE(Cish2)
                  # P-value, Cis-SNPs
                  line = c(substring(file_name, 1, nchar(file_name) - 4), 
                           hsq[1,2], hsq[2,2], hsq[3,2], hsq[4,2],
                           hsq[1,3], hsq[2,3], hsq[3,3], hsq[4,3],
                           hsq_log[5,2], hsq_log[7,2]
                  )
                  return(line)
            }))
colnames(hsqs) = c("ID", "V_Cis", "V_Res", "V_Phe", "Cis_h2",
                   "SE_V_Cis", "SE_V_Res", "SE_V_Phe", "SE_Cis_h2",
                   "P_val", "Cis_SNPs"
)
rownames(hsqs) = hsqs[,"ID"]
hsqs = as.data.frame(hsqs)
for (i in 2:10) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/all_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)
