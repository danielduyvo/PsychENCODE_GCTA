#!/usr/bin/env Rscript

PROJECT="EUR_SPC_HRC_bksk_0.025"

# Combine the HSQ files and return statistics
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/hsqs", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                  hsq = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=TRUE, nrows=7, fill=TRUE)
                  hsq[7,] = c("V(G)/Vp", hsq[8,1], hsq[8,2])
                  hsq = hsq[-8,]
                  hsq_log = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=FALSE, skip=9, fill=TRUE)
                  # Gene ID,
                  # Variance from bK GRM, Variance from SNP, Residual Variance, 
                  # Phenotypic Variance, 
                  # narrowh2 - SNPh2, SNPh2, narrowh2
                  # SE(Variance from bK GRM), SE(Variance from SNP), SE(Residual Variance),
                  # SE(Phenotypic Variance),
                  # SE(narrowh2 - SNPh2), SE(SNPh2), SE(narrowh2),
                  # P_val
                  line = c(substring(file_name, 1, nchar(file_name) - 4), 
                           hsq[1,2], hsq[2,2], hsq[3,2], hsq[4,2], 
                           hsq[5,2], hsq[6,2], hsq[7,2],
                           hsq[1,3], hsq[2,3], hsq[3,3], hsq[4,3],
                           hsq[5,3], hsq[6,3], hsq[7,3],
                           hsq_log[5,2]
                  )
                  return(line)
            }))
colnames(hsqs) = c("ID", "V_bK", "V_SNP", "Res", "Phe", "bKh2", "SNPh2", "narrowh2",
                   "SE_V_bk", "SE_V_SNP", "SE_Res", "SE_Phe", "SE_bKh2", "SE_SNPh2", "SE_narrowh2",
                   "P_val"
)
rownames(hsqs) = hsqs[,"ID"]
hsqs = as.data.frame(hsqs)
for (i in 2:16) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/all_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

avg_snp_h2 = mean(hsqs$SNPh2)
narrow_sense_h2 = mean(hsqs$narrowh2)
