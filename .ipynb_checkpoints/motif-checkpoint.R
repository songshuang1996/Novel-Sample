#library(MotifDb)
library(universalmotif)
library(ggplot2)
##读取共用数据库中lexA的结果
#motifs <- filter_motifs(MotifDb, organism = c("Athaliana", "Celegans"))



# motif<-read_meme('~/Desktop/meme.txt')
# motif1 <- filter_motifs(motif, icscore = 0)
# dedup_motif<-merge_similar(motif1)
# logplot<-view_motifs(dedup_motif)
# #ggsave('~/Desktop/DR_motifs.pdf',plot=logplot,device = 'pdf',width = 250, height = 2000, limitsize = FALSE,units = "cm")
# newmotif=convert_type(motif1,'PPM')
# newotif <- newmotif["motif"]

print(getwd())
#循环读入所有的motif
rm(list=ls())
motiftogether=c()
#setwd('/data/songshuang/server_find_novel_regulator/99_User_project/2023/3_meme_result')
organs=list.files()

for (organ in organs) {
  #读取motif
  motif=read_meme(sprintf('%s/meme.txt',organ))
  if (length(motif)==0) {
    next
  }
  if (length(motif)==1) {
    motif=c(motif)
  }
  #重命名
  for (number in 1:length(motif)) {
    motif[[number]]['name']=sprintf('pal-%s-%s-%s',organ,length(motif),number)
  }
  motiftogether=c(motiftogether,motif)
}

motiftogether=filter_motifs(motiftogether,width = 10,icscore = 15)
motiftogether=trim_motifs(motiftogether,min.ic=0.8)
motiftogether=filter_motifs(motiftogether,width = 10,icscore = 15)

motiftogether<-merge_similar(motiftogether,threshold = 0.6,nthreads = 32)
logplot<-view_motifs(motiftogether)
ggsave('novel_motif.png',plot=logplot,width = 50, height = length(motiftogether)*8, units = "cm",device = 'png',limitsize = FALSE)





