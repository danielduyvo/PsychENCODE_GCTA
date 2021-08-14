#!/usr/bin/env Rscript

PROJECT="EUR_SPC_HRC_whole"

# Combine the HSQ files and return statistics
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/hsqs", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                  hsq = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=TRUE, nrows=4, fill=TRUE)
                  hsq_log = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=FALSE, skip=5, fill=TRUE)
                  # Gene ID, 
                  # Variance from SNP, Residual Variance, Phenotypic Variance, SNPh2,
                  # SE(Variance from SNP), SE(Residual Variance), SE(Phenotypic Variance), SE(SNPh2)
                  # P-value
                  line = c(substring(file_name, 1, nchar(file_name) - 4), 
                           hsq[1,2], hsq[2,2], hsq[3,2], hsq[4,2],
                           hsq[1,3], hsq[2,3], hsq[3,3], hsq[4,3],
                           hsq_log[5,2]
                  )
                  return(line)
            }))
colnames(hsqs) = c("ID", "V_SNP", "V_Res", "V_Phe", "SNP_h2",
                   "SE_V_SNP", "SE_V_Res", "SE_V_Phe", "SE_SNP_h2",
                   "P_val"
)
rownames(hsqs) = hsqs[,"ID"]
hsqs = as.data.frame(hsqs)
for (i in 2:10) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/all_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

avg_snp_h2 = mean(hsqs[, "SNPh2"]) # 0.06808974
