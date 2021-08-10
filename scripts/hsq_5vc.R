#!/usr/bin/env Rscript
PROJECT="AFR_SPC_HRC_bksk_1Mbase"
BKSK="AFR_SPC_HRC_bksk_0.025"

# Grabbing data from HSQs produced by genes with SNPs in the cis window
temp = unlist(list.files(path=paste("data/", PROJECT, "/output/hsqs", sep=""), pattern="*.hsq"))
hsqs = t(sapply(X=temp, function(file_name) {
                    hsq = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=TRUE, nrows=12, fill=TRUE)
                    hsq[11,] = c("V(G)/Vp", hsq[12,1], hsq[12,2])
                    hsq = hsq[-12,]
                    hsq_log = read.table(paste('data/', PROJECT, '/output/hsqs/', file_name, sep=''), header=FALSE, skip=13, fill=TRUE)
                    # Gene ID,
                    # Variance from Cis bK, Variance from Cis SNP,
                    # Variance from Trans bK, Variance from Trans SNP,
                    # Residual Variance, Phenotypic Variance, 
                    # cis bK h2, cis SNP h2
                    # trans bK h2, trans SNP h2,
                    # narrow sense h2
                    # SE(Variance from Cis bK), SE(Variance from Cis SNP),
                    # SE(Variance from Trans), SE(Variance from Trans SNP),
                    # SE(Residual Variance), SE(Phenotypic Variance),
                    # SE(cis bK h2), SE(cis SNP h2),
                    # SE(transh2), SE(trans SNP h2),
                    # SE(narrow sense h2)
                    # P-value
                    line = c(substring(file_name, 1, nchar(file_name) - 4), 
                             hsq[1,2], hsq[2,2], 
                             hsq[3,2], hsq[4,2], 
                             hsq[5,2], hsq[6,2], 
                             hsq[7,2], hsq[8,2],
                             hsq[9,2], hsq[10,2],
                             hsq[11,2],
                             hsq[1,3], hsq[2,3], 
                             hsq[3,3], hsq[4,3], 
                             hsq[5,3], hsq[6,3], 
                             hsq[7,3], hsq[8,3],
                             hsq[9,3], hsq[10,3],
                             hsq[11,3],
                             hsq_log[5,2]
                    )
                    return(line)
}))
colnames(hsqs) = c("ID", 
                   "V_Cis_bK", "V_Cis_SNP", 
                   "V_Trans_bK", "V_Trans_SNP", 
                   "V_Res", "V_Phe",
                   "Cis_bk_h2", "Cis_SNP_h2",
                   "Trans_bk_h2", "Trans_SNP_h2",
                   "Narrow_h2",
                   "SE_V_Cis_bK", "SE_V_Cis_SNP",
                   "SE_V_Trans_bK", "SE_V_Trans_SNP", 
                   "SE_V_Res", "SE_V_Phe",
                   "SE_Cis_bK_h2", "SE_Cis_SNP_h2",
                   "SE_Trans_bK_h2", "SE_Trans_SNP_h2",
                   "SE_Narrow_h2",
                   "P_val"
)
rownames(hsqs) = hsqs[,"ID"]
hsqs = as.data.frame(hsqs)
for (i in 2:24) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/all_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

system(paste("export PROJECT=", PROJECT, "; ./scripts/missing_hsq_1Mbase_helper.sh", sep=""))

# Process HSQs from genes missing SNPs in the cis windows
temp = read.table(paste("data/", PROJECT, "/output/results/missing.txt", sep=""), header=FALSE)$V1
hsqs = t(sapply(X=temp, function(file_name) {
                    hsq = read.table(paste('data/', BKSK, '/output/hsqs/', file_name, sep=''), header=TRUE, nrows=7, fill=TRUE)
                    hsq[7,] = c("V(G)/Vp", hsq[8,1], hsq[8,2])
                    hsq = hsq[-8,]
                    hsq_log = read.table(paste('data/', BKSK, '/output/hsqs/', file_name, sep=''), header=FALSE, skip=9, fill=TRUE)

                    # Gene ID,
                    # Variance from Cis bK, Variance from Cis SNP,
                    # Variance from Trans bK, Variance from Trans SNP,
                    # Residual Variance, Phenotypic Variance, 
                    # cis bK h2, cis SNP h2
                    # trans bK h2, trans SNP h2,
                    # narrow sense h2
                    # SE(Variance from Cis bK), SE(Variance from Cis SNP),
                    # SE(Variance from Trans), SE(Variance from Trans SNP),
                    # SE(Residual Variance), SE(Phenotypic Variance),
                    # SE(cis bK h2), SE(cis SNP h2),
                    # SE(transh2), SE(trans SNP h2),
                    # SE(narrow sense h2)
                    # P-value
                    line = c(substring(file_name, 1, nchar(file_name) - 4), 
                             0, 0,
                             hsq[1,2], hsq[2,2], 
                             hsq[3,2], hsq[4,2], 
                             0, 0,
                             hsq[5,2], hsq[6,2],
                             hsq[7,2],
                             0, 0,
                             hsq[1,3], hsq[2,3], 
                             hsq[3,3], hsq[4,3], 
                             0, 0,
                             hsq[5,3], hsq[6,3],
                             hsq[7,3],
                             hsq_log[5,2]
                    )
                    return(line)
            }))
colnames(hsqs) = c("ID", 
                   "V_Cis_bK", "V_Cis_SNP", 
                   "V_Trans_bK", "V_Trans_SNP", 
                   "V_Res", "V_Phe",
                   "Cis_bk_h2", "Cis_SNP_h2",
                   "Trans_bk_h2", "Trans_SNP_h2",
                   "Narrow_h2",
                   "SE_V_Cis_bK", "SE_V_Cis_SNP",
                   "SE_V_Trans_bK", "SE_V_Trans_SNP", 
                   "SE_V_Res", "SE_V_Phe",
                   "SE_Cis_bK_h2", "SE_Cis_SNP_h2",
                   "SE_Trans_bK_h2", "SE_Trans_SNP_h2",
                   "SE_Narrow_h2",
                   "P_val"
)
rownames(hsqs) = hsqs[,"ID"]
hsqs = as.data.frame(hsqs)
for (i in 2:24) {
    hsqs[,i] = as.numeric(hsqs[,i])
}
write.table(hsqs, paste("data/", PROJECT, "/output/results/no_cis_variance.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

# Combine HSQs from both sets to generate statistics for the whole dataset
hsqs = read.table(paste('data/', PROJECT, '/output/results/all_variance.txt', sep=''), header=TRUE)
no_cis_hsqs = read.table(paste('data/', PROJECT, '/output/results/no_cis_variance.txt', sep=''), header=TRUE)
all_hsqs = rbind(hsqs, no_cis_hsqs)

