library(magrittr)

setwd('C:\\Users\\ftroo\\Desktop\\Lab\\Research\\5_SCseq_LCM\\20210427_Generate_data_for_MetaCell')

RDSpath = 'C:\\Users\\ftroo\\Desktop\\Lab\\Research\\5_SCseq_LCM\\20210412_Reorganize_the_data\\'

Ptr_SCseq_UMI_matrix = readRDS(paste0(RDSpath,'RDS_Ptr_SCseq_UMI_matrix.rds'))
Lch_SCseq_UMI_matrix = readRDS(paste0(RDSpath,'RDS_Lch_SCseq_UMI_matrix.rds'))

Ptr_SCseq_UMI_table = as.data.frame(Ptr_SCseq_UMI_matrix)
write.table(Ptr_SCseq_UMI_table,'PoT_SCseq_UMI_table.tsv',sep='\t',quote=F)
cat(rownames(Ptr_SCseq_UMI_table),file='PoT_CDS.txt',sep='\n')

Lch_SCseq_UMI_table = as.data.frame(Lch_SCseq_UMI_matrix)
write.table(Lch_SCseq_UMI_table,'LiC_SCseq_UMI_table.tsv',sep='\t',quote=F)
cat(rownames(Lch_SCseq_UMI_table),file='LiC_CDS.txt',sep='\n')

