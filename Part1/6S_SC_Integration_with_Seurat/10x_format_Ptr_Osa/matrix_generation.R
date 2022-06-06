library(magrittr)
library(dplyr)
library(Matrix)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210527_10x_format_Ptr_Osa')


#Read in the gff table for CDS2gene ID conversion
Osa_gff_table = read.table('Oryza_sativa.IRGSP-1.0.51.chr.gff3',
                           sep='\t',header=F,skip=20,quote='')
dim(Osa_gff_table) # 545322 9
Osa_CDS_info = filter(Osa_gff_table,V3=='CDS')$V9 %>% unique
length(Osa_CDS_info) # 42355
Osa_CDS_ID = Osa_CDS_info %>% strsplit(';') %>% sapply(extract,1) %>% sub('ID=CDS:','',.)
Osa_transcript_ID = Osa_CDS_info %>% strsplit(';') %>% sapply(extract,2) %>% sub('Parent=transcript:','',.)
all(Osa_CDS_ID==Osa_transcript_ID) # FALSE: not all CDS_ID==transcript_ID
Osa_CDS2transcriptID = data.frame(cbind(Osa_CDS_ID,Osa_transcript_ID))

Osa_mRNA_info = filter(Osa_gff_table,V3=='mRNA')$V9 %>% unique
length(Osa_mRNA_info) # 44754
Osa_transcript_ID = Osa_mRNA_info %>% strsplit(';') %>% sapply(extract,1) %>% sub('ID=transcript:','',.)
Osa_gene_ID = Osa_mRNA_info %>% strsplit(';') %>% sapply(extract,2) %>% sub('Parent=gene:','',.)
Osa_transcript2geneID = data.frame(cbind(Osa_transcript_ID,Osa_gene_ID))

all(Osa_CDS2transcriptID$Osa_transcript_ID %in% Osa_transcript2geneID$Osa_transcript_ID)
Osa_CDS2geneID = full_join(Osa_CDS2transcriptID,Osa_transcript2geneID,by='Osa_transcript_ID')
dim(Osa_CDS2geneID) # 44754 3

##We only keep one CDS ID for each gene ID
##There are two types of CDS ID: "cds-[BAC|CCA]XXXXX.X" and "OsXXtXXXXXXX-XX"
Osa_CDS2geneID$Osa_CDS_ID %>% substr(1,5) %>% table # 109 with prefix "cds"

part1 = filter(na.omit(Osa_CDS2geneID),grepl('cds',Osa_CDS_ID))
part2 = filter(na.omit(Osa_CDS2geneID),!grepl('cds',Osa_CDS_ID))

(length(unique(part1$Osa_gene_ID))==109) #TRUE: for those CDS prefix as "cds", each gene have single isoform

#Handle with those with prefix as "OsXXt"
all(sapply(part2$Osa_CDS_ID,nchar) == 15) # TRUE
all(substr(part2$Osa_CDS_ID,1,12) == sub('g','t',part2$Osa_gene_ID)) # TRUE

##Keep the CDS with the smallest isoform id (CDS_subID)
Osa_unqiue_gene_ID = unique(part2$Osa_gene_ID)
Osa_CDS_subID = substr(part2$Osa_CDS_ID,14,15)
table(Osa_CDS_subID) # Max:09
selected_CDS_ID = c()
selected_gene_ID = c()
for(i in 0:9){
  candidate_gene_ID = part2$Osa_gene_ID[Osa_CDS_subID==sprintf('%02d',i)]
  new_gene_ID = setdiff(candidate_gene_ID,selected_gene_ID)
  if(length(new_gene_ID)>0){
    selected_gene_ID %<>% c(new_gene_ID)
    selected_CDS_ID %<>% c(paste0(sub('g','t',new_gene_ID),sprintf('-%02d',i)))
  }
  if(all(Osa_unqiue_gene_ID %in% selected_gene_ID)){
    print(i)
    break
  }
}

total_selected_CDS_ID = c(part1$Osa_CDS_ID,selected_CDS_ID)
Osa_selected_CDS2geneID = filter(Osa_CDS2geneID,Osa_CDS_ID%in%total_selected_CDS_ID)
dim(Osa_selected_CDS2geneID) # 35775 3


#Prepare the ortholog cluster table
ortho_cluster_table = readRDS('RDS_ortho_cluster_table.rds')

Ptr_ortho_cluster_table = filter(ortho_cluster_table,speciesName=='PoT')

Osa_ortho_cluster_table = filter(ortho_cluster_table,
                                 speciesName=='OrS',
                                 geneName%in%Osa_selected_CDS2geneID$Osa_CDS_ID)
Osa_ortho_cluster_table$geneName = with(Osa_selected_CDS2geneID,
                                        Osa_gene_ID[match(Osa_ortho_cluster_table$geneName,Osa_CDS_ID)])
dim(Osa_ortho_cluster_table) # 35775 3


#Read the UMI matrix
##Ptr
Ptr_SCseq_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')

##Osa
read10xMatrix = function(Path='2021_Zhang\\osRoot1\\filtered_feature_bc_matrix\\'){
  mat = readMM(paste0(Path,'matrix.mtx.gz'))
  barcode.names = read.delim(paste0(Path,'barcodes.tsv.gz'),header=F,stringsAsFactors=F)
  feature.names = read.delim(paste0(Path,'features.tsv.gz'),header=F,stringsAsFactors=F)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  return(mat)
}

Osa_2021R1_SCseq_UMI_matrix = read10xMatrix('2021_Zhang\\osRoot1\\filtered_feature_bc_matrix\\')
Osa_2021R2_SCseq_UMI_matrix = read10xMatrix('2021_Zhang\\osRoot2\\filtered_feature_bc_matrix\\')
dim(Osa_2021R1_SCseq_UMI_matrix) # 38978 15296
dim(Osa_2021R2_SCseq_UMI_matrix) # 38978 29246


unwanted_geneID = filter(Osa_gff_table,V3 %in% c('ncRNA_gene','pseudogene'))$V9 %>% 
  strsplit(';') %>% sapply(extract,1) %>% sub('ID=gene:','',.)

noCDS_geneID = filter(Osa_CDS2geneID,is.na(Osa_CDS_ID))$Osa_gene_ID %>% unique

setequal(rownames(Osa_2021R1_SCseq_UMI_matrix),
         c(Osa_selected_CDS2geneID$Osa_gene_ID,
           unwanted_geneID,noCDS_geneID)) # TRUE


Osa_2021R1_SCseq_UMI_matrix %<>% extract(Osa_selected_CDS2geneID$Osa_gene_ID,)
Osa_2021R2_SCseq_UMI_matrix %<>% extract(Osa_selected_CDS2geneID$Osa_gene_ID,)
dim(Osa_2021R1_SCseq_UMI_matrix) # 35775 15296
dim(Osa_2021R2_SCseq_UMI_matrix) # 35775 29246

saveRDS(Osa_2021R1_SCseq_UMI_matrix,
        'RDS_Osa_2021R1_SCseq_UMI_matrix.rds')
saveRDS(Osa_2021R2_SCseq_UMI_matrix,
        'RDS_Osa_2021R2_SCseq_UMI_matrix.rds')






#Calculate the ortholog cluster expression
overlap_cluster = intersect(unique(Ptr_ortho_cluster_table$clusterID),
                            unique(Osa_ortho_cluster_table$clusterID))
length(overlap_cluster) # 9800

calculate_ortho_UMI_matrix = function(SCseq_UMI_matrix = Ptr_SCseq_UMI_matrix,
                                      ortho_cluster_table = Ptr_ortho_cluster_table,
                                      cal_cluster = overlap_cluster,
                                      output_prefix='Ptr_',
                                      save_RDS=F){
  geneID_list_of_each_cluster = sapply(cal_cluster,
                                       function(cc) with(ortho_cluster_table,
                                                         geneName[clusterID==cc]))
  matrix_geneID = rownames(SCseq_UMI_matrix)
  cluster_indicator_matrix = Matrix(0,
                                    nrow=length(cal_cluster),
                                    ncol=length(matrix_geneID),
                                    dimnames=list(paste0('Cluster_',cal_cluster),
                                                  matrix_geneID))
  for(i in seq_along(geneID_list_of_each_cluster)){
    selected_genes = geneID_list_of_each_cluster[[i]]
    cluster_indicator_matrix[i,selected_genes] = 1
  }
  SCseq_ortho_cluster_UMI_matrix = as.matrix(cluster_indicator_matrix%*%Matrix(SCseq_UMI_matrix))
  if(save_RDS) saveRDS(SCseq_ortho_cluster_UMI_matrix,
                       paste0('RDS_',output_prefix,'SCseq_ortho_cluster_UMI_matrix.rds'))
  return(SCseq_ortho_cluster_UMI_matrix)
}


##Ptr
Ptr_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Ptr_SCseq_UMI_matrix,
                             ortho_cluster_table = Ptr_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Ptr_',
                             save_RDS=F)
dim(Ptr_SCseq_ortho_cluster_UMI_matrix) # 9800 4705


##Osa: 2021_Zhang
Osa_2021R1_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Osa_2021R1_SCseq_UMI_matrix,
                             ortho_cluster_table = Osa_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Osa_2021R1_',
                             save_RDS=F)
dim(Osa_2021R1_SCseq_ortho_cluster_UMI_matrix) # 9800 15296
Osa_2021R2_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Osa_2021R2_SCseq_UMI_matrix,
                             ortho_cluster_table = Osa_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Osa_2021R2_',
                             save_RDS=F)
dim(Osa_2021R2_SCseq_ortho_cluster_UMI_matrix) # 9800 29246



#Prepare the 10x format files for Seurat
output_as_10x_format = function(matrix,prefix){
  
  gz_temp = gzfile(paste0(prefix,'features.tsv.gz'),'w')
  foo = rownames(matrix)
  write.table(data.frame(foo,foo,'Gene Expression'),
              gz_temp,
              sep='\t',row.names=F,col.names=F,quote=F)
  close(gz_temp)
  
  gz_temp = gzfile(paste0(prefix,'barcodes.tsv.gz'),'w')
  write.table(data.frame(paste0(prefix,colnames(matrix))),
              gz_temp,
              sep='\t',row.names=F,col.names=F,quote=F)
  close(gz_temp)
  
  writeMM(Matrix(matrix),paste0(prefix,'matrix.mtx'))
  system(paste0('gzip ',paste0(prefix,'matrix.mtx')))
  
  return(0)
}


# Ptr_SCseq_ortho_cluster_UMI_matrix = readRDS('RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')

output_as_10x_format(Ptr_SCseq_ortho_cluster_UMI_matrix,'Ptr_')
output_as_10x_format(Osa_2021R1_SCseq_ortho_cluster_UMI_matrix,'Osa_2021R1_')
output_as_10x_format(Osa_2021R2_SCseq_ortho_cluster_UMI_matrix,'Osa_2021R2_')


