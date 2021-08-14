#!/usr/bin/env Rscript
PROJECT="EUR_SPC_HRC_chrs"

# Grabbing data from HSQs
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/hsqs", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                    hsq = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=TRUE, nrows=7, fill=TRUE)
                    hsq[7,] = c("V(G)/Vp", hsq[8,1], hsq[8,2])
                    hsq = hsq[-8,]
                    hsq_log = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=FALSE, skip=9, fill=TRUE)
                    # Gene ID,
                    # Variance from Cis, Variance from Trans, Residual Variance, 
                    # Phenotypic Variance, 
                    # cish2, transh2, SNPh2
                    # SE(Variance from Cis), SE(Variance from Trans), SE(Residual Variance),
                    # SE(Phenotypic Variance),
                    # SE(cish2), SE(transh2), SE(SNPh2)
                    # P-value
                    line = c(substring(file_name, 1, nchar(file_name) - 4), 
                             hsq[1,2], hsq[2,2], hsq[3,2], hsq[4,2], 
                             hsq[5,2], hsq[6,2], hsq[7,2],
                             hsq[1,3], hsq[2,3], hsq[3,3], hsq[4,3],
                             hsq[5,3], hsq[6,3], hsq[7,3],
                             hsq_log[5,2]
                    )
                    return(line)
}))
colnames(hsqs) = c("ID", "V_Cis", "V_Trans", "V_Res", "V_Phe",
                   "Cis_h2", "Trans_h2", "SNP_h2",
                   "SE_V_Cis", "SE_V_Trans", "SE_V_Res", "SE_V_Phe",
                   "SE_Cis_h2", "SE_Trans_h2", "SE_SNP_h2",
                   "P_val"
)
rownames(hsqs) = hsqs[,"ID"]
hsqs = as.data.frame(hsqs)
for (i in 2:16) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/all_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)
