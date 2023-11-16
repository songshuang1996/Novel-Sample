library(universalmotif)
library(ggplot2)
library(stringr)

#read all motifs
rm(list=ls())
motiftogether=c()
organs=list.files()

for (organ in organs) {
  #read in all meme results

  motif=read_meme(sprintf('%s/meme.txt',organ))
  
  if (length(motif)==0) {
    next
  }
  if (length(motif)==1) {
    motif=c(motif)
  }
  #rename
  for (number in 1:length(motif)) {
    motif[[number]]['name']=sprintf('%s-%s-%s',organ,length(motif),number)
  }
  motiftogether=c(motiftogether,motif)
}


#cleaning the motifs 
motiftogether=filter_motifs(motiftogether,width = 14,icscore = 15)
motiftogether=trim_motifs(motiftogether,min.ic=0.6)
motiftogether=filter_motifs(motiftogether,width = 14,icscore = 15)
#merge similar motifs
motiftogether1<-merge_similar(motiftogether,threshold = 0.6,nthreads = 64)
motiftogether1=trim_motifs(motiftogether1,min.ic=0.6)

count_motif=0
gene_count=c()
motiftogether2=motiftogether1

#give a summary of the results
if (length(motiftogether1)>1){
  for (single_motif in motiftogether1) {
    count_motif=count_motif+1
    motif_name=motiftogether1[[count_motif]]['name']
    protein_id_all=str_replace_all(string = motif_name,pattern = "-\\d-\\d",replacement = "")
    occur_times=str_count(motif_name,'/')[1]+1
    split_data=str_split(protein_id_all,'/|-')
    
    motiftogether2[[count_motif]]['name']=sprintf('Motif %d - Frequency of occurrence:%d',count_motif,occur_times)
    #how many time does a gene in combinations
    count_data=sort(table(split_data),decreasing = T)
    write.table(count_data,sprintf('Motif_%d',count_motif),row.names = F,col.names = F,sep = ',')
}
}else{
  motif_name=motiftogether1[[1]]['name']
  protein_id_all=str_replace_all(string = motif_name,pattern = "-\\d-\\d",replacement = "")
  occur_times=str_count(motif_name,'/')[1]
  split_data=str_split(protein_id_all,'/|-')
  
  motiftogether2[[1]]['name']=sprintf('Motif %d - Frequency of occurrence:%d',1,occur_times)
  #how many time does a gene in combinations
  count_data=sort(table(split_data),decreasing = T)
  write.table(count_data,sprintf('Motif_1'),row.names = F,col.names = F,sep = ',')
}


logplot<-view_motifs(motiftogether2)
ggsave('novel_motif.png',plot=logplot,width = 50, height = length(motiftogether2)*8, units = "cm",device = 'png',limitsize = FALSE)





