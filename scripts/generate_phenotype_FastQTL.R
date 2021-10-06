#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

bed_filename = args[1]
population = args[2]
pheno = args[3]

bed_file = read.table(pipe(paste('zcat', bed_filename, '| sed -e 1s/#//', sep=' ')), 
                      header = TRUE, sep='\t')
bed_file[, 'Chr'] = sapply(X = bed_file[, 'Chr'], FUN = function(chr_str) {
                           return(gsub('[^0-9]', '', chr_str))
                      })
bed_file[bed_file[, 'Chr'] == '', 'Chr'] = '0'
fam_id = rep(x = 0, times=dim(bed_file)[2] - 4)
sample_id = scan(text = readLines(bed_filename, 1), what = "", 
                 quiet = TRUE)[-(1:4)]
phenotype_id = bed_file[,'ID']
phenotypes = bed_file[,5:dim(bed_file)[2]]
phenotypes = t(as.matrix(phenotypes))
colnames(phenotypes) = phenotype_id
rownames(phenotypes) = sample_id

phenotype_info = cbind(phenotype_id, bed_file[,c('Chr', 'start', 'end')])

write.table(phenotype_info, paste('fin_RNA_seq/', population, '_', pheno, '/phenotype_ids', sep=''), 
            col.names=FALSE, row.names=FALSE, quote=FALSE)

pheno_file = cbind(fam_id, sample_id, phenotypes)
write.table(pheno_file, paste('fin_RNA_seq/', population, '_', pheno, '/phenotype', sep=''), 
            sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
