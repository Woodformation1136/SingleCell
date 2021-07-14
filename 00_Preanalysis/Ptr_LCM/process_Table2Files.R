library(magrittr)

Args = commandArgs(trailingOnly=T)

inTable = read.csv(Args[1])

dir.create('./outputSampleList', showWarnings=F)
for(i in seq(nrow(inTable))){
  Out = data.frame(R1=inTable[i,'Read_1'],
                   R2=inTable[i,'Read_2'],
                   ID=paste0('ID:',inTable[i,'Group']))
  write.table(Out,file=paste0('./outputSampleList/',inTable[i,'Group']),
              col.names=F,row.names=F,quote=F,sep='\t')
}


