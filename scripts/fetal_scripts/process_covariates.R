library(data.table)
cov = read.table('/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/50HCP_cov_trimester1.txt',
                 header=FALSE,
                 sep='\t'
)

t_cov = transpose(cov[,-1]) # Omit covariate name

t_cov = cbind(rep(x=0, times=dim(t_cov)[1]), t_cov)
colnames(t_cov) = c('fam_id', cov[,1])

# Generate covariates
qual_cov=t_cov[,c("fam_id", "id", "sex")]

# Generate qcovariates
quant_cov=t_cov[,-8]

write.table(qual_cov, 'fin_RNA_seq/Trimester_1_EUR_50_HPC_gene/qual_covars.txt',
            quote=FALSE,
            sep='\t',
            col.names=FALSE,
            row.names=FALSE
)
write.table(quant_cov, 'fin_RNA_seq/Trimester_1_EUR_50_HPC_gene/quant_covars.txt',
            quote=FALSE,
            sep='\t',
            col.names=FALSE,
            row.names=FALSE
)

cov = read.table('/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/50HCP_cov_trimester2.txt',
                 header=FALSE,
                 sep='\t'
)

t_cov = transpose(cov[,-1]) # Omit covariate name

t_cov = cbind(rep(x=0, times=dim(t_cov)[1]), t_cov)
colnames(t_cov) = c('fam_id', cov[,1])

# Generate covariates
qual_cov=t_cov[,c("fam_id", "id", "sex")]

# Generate qcovariates
quant_cov=t_cov[,-8]

write.table(qual_cov, 'fin_RNA_seq/Trimester_2_EUR_50_HPC_gene/qual_covars.txt',
            quote=FALSE,
            sep='\t',
            col.names=FALSE,
            row.names=FALSE
)
write.table(quant_cov, 'fin_RNA_seq/Trimester_2_EUR_50_HPC_gene/quant_covars.txt',
            quote=FALSE,
            sep='\t',
            col.names=FALSE,
            row.names=FALSE
)
