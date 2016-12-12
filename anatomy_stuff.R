
#load the registration results
load('~/GitHub/St-dat/Registration_results/160919_SE2_St&Th_ID3_E1_H&E_AJu.RData')

#displays what you have in workingspace
ls()

#features = spots of Cy3
#cell bodies = cell bodies
#regi registration object


#install here http://www.wholebrainsoftware.org/cms/installing-wholebrain-on-mac-osx/
library(wholebrain)

#wholebrain package
dataset<-inspect.registration(regi,features, forward.warps=TRUE)
regi1<-regi
head(dataset)

load('~/GitHub/St-dat/Registration_results/160919_SE2_St&Th_ID3_E2_H&E_AJu.RData')
dataset2<-inspect.registration(regi,features, forward.warps=TRUE)
regi2<-regi


#default R
sort(table(dataset$name))

sort(table(dataset$acronym))


plot(dataset$x, dataset$y, ylim=c(regi$transformationgrid$height,0), xlim=c(0,regi$transformationgrid$width), col=0 , asp=1)

polygon(rep(c(0,regi1$transformationgrid$width), each=2),c(c(regi1$transformationgrid$height,0), rev(c(regi1$transformationgrid$height,0))), col='black')

points(dataset$x, dataset$y, pch=21, bg='orange')
points(dataset2$x, dataset2$y, pch=21, bg='purple')


plot(dataset$ML, dataset$DV, col=0 , asp=1)

index1<-which(!is.na(dataset$acronym))
index2<-which(!is.na(dataset2$acronym))

points(dataset$ML[index1], dataset$DV[index1], pch=21, bg='orange')
points(dataset2$ML[index2], dataset2$DV[index2], pch=21, bg='purple')


read.ST<-function(filepath='adjusted_stdata.tsv'){
dat<-read.table(filepath, sep='\t', header=TRUE, row.names=1)

xcoord<-unlist(strsplit(row.names(dat),'x'))
  ycoord<-xcoord
  
xcoord<-as.numeric(as.vector(xcoord[seq(1, length(xcoord), by=2)]))

ycoord<-as.numeric(as.vector(ycoord[seq(2, length(ycoord), by=2)]))


dat<-data.frame(y=ycoord, x=xcoord, dat)

data<-dat

return(data)
}

plot(as.numeric(dat1[,1]), as.numeric(dat1[,2]))

dat1<-read.ST('~/GitHub/St-dat/Sequencing_reads/E1_ID3_S3_stdata.tsv')
dat1.adj<-read.ST('~/GitHub/St-dat/Sequencing_reads/E2_ID3_S3_stdata_adj.tsv')

dat2<-read.ST('~/GitHub/St-dat/Sequencing_reads/E2_ID3_S3_stdata.tsv')

dataset<-dataset[-c(which.min(dataset$y), which.max(dataset$y)),]

quartz()
plot(as.numeric(as.vector(dat1[,1])), as.numeric(as.vector(dat1[,2])), cex=3, ylim=c(32, 2), xlim=c(34,2))
text(as.numeric(as.vector(dat1[,1])), as.numeric(as.vector(dat1[,2])), labels=row.names(dat1), cex=0.5)
points(as.numeric(as.vector(dat1.adj[,1])), as.numeric(as.vector(dat1.adj[,2])), cex=3, pch=21, bg='red')

plot(dataset$x, dataset$y, bg=dataset$color, pch=21)


names(dat1)[3:ncol(dat1)]%in%

