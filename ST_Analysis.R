
#load file
read.ST<-function(filepath='adjusted_stdata.tsv'){
dat<-read.table(filepath, sep='\t', header=TRUE, row.names=1)

xcoord<-unlist(strsplit(row.names(dat),'x'))
  ycoord<-xcoord
  
xcoord<-xcoord[seq(1, length(xcoord), by=2)]

ycoord<-ycoord[seq(2, length(ycoord), by=2)]


dat<-data.frame(y=ycoord, x=xcoord, dat)

data<-dat

return(data)
}


dat1<-read.ST('/Users/annyschmotzke/Documents/Work/ST-data/Layer_12_rep1/Data/adjusted_stdata.tsv')
dat2<-read.ST('/Users/annyschmotzke/Documents/Work/ST-data/Layer_12_rep2/Data/adjusted_stdata.tsv')

data<-rbind(dat1, dat2)
data<-data.frame(replicate=c(rep(1,nrow(dat1)), rep(2,nrow(dat2))), data)

plot(ycoord, xcoord, ylim=c(35,0), xlim=c(0,33))

dat<-log10(dat+1)
dispersion<- apply(dat, 2, var)
plot(apply(dat, 2, mean), dispersion)
abline(h=quantile(dispersion, 0.95))
[,which(dispersion>quantile(dispersion, 0.99))]


d <- dist(data, method = "manhattan") # distance matrix
fit <- hclust(d, method="ward") 

numbers.clusters<-15

groups <- cutree(fit, k=15) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 


par(mfrow=c(1,2))
cluster.palette<-rgb(runif(100),runif(100),runif(100))

plot(ycoord, xcoord, ylim=c(35,0), xlim=c(0,33), pch=16, col=cluster.palette[groups])

plot(fit, col='gray') # display dendogram

rect.hclust(fit, k=numbers.clusters, border=cluster.palette[unique(groups)])

