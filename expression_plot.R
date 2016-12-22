dataset$ID<-paste(dataset$array.y, dataset$array.x, sep='x')

dataset2$ID<-paste(dataset2$array.y, dataset2$array.x, sep='x')

dataset3<-dataset2
dataset1<-dataset

dataset3<-dataset3[which(dataset3$ID%in%data2$ID),]
dataset3<-dataset3[order(dataset3$ID),]

dataset1 <-dataset1[which(dataset1$ID%in%data1$ID),]
dataset1<-dataset1[order(dataset$ID),]

data1 <-data1[order(data1$ID),]
data2 <-data2[order(data2$ID),]


plot.gene()
plot.gene(ensemb ='ENSMUSG00000020178.5', sym='Adora2a')
plot.gene(ensemb = 'ENSMUSG00000032532.7', sym='Cck')
plot.gene(ensemb='ENSMUSG00000061762.12', sym='Tac1')
plot.gene(ensemb='ENSMUSG00000004151.16', sym='Etv1')
plot.gene(ensemb='ENSMUSG00000024907.6', sym='Gal')
plot.gene(ensemb='ENSMUSG00000092341.2', sym='Malat1')

plot.gene<-function(ensemb ='ENSMUSG00000021478.5', sym='Drd1'){
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = 1)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
plot.new()
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")

rect(0, levels[-length(levels)], 1, levels[-1L],col=rev(col),border=rev(col)) 
axis(4)
box()

mar <- mar.orig
   mar[4L] <- 1
   par(mar = mar)


plot(dataset$ML, dataset$DV, col=0 , asp=1, xlim=c(-5,0), ylab='Dorso-ventral (mm)', xlab='Medio-lateral (mm)', axes=F, main=paste(sym, '\n', ensemb))
abline(h=c(3:-10), col='lightblue')
abline(v=c(8:-8), col='lightblue')

scale.factor<-mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height,registration$transformationgrid$width) )
 numPaths<-registration$atlas$numRegions
 outlines<-registration$atlas$outlines


lapply(1:numPaths, function(x){polygon(stereotactic.coordinates(outlines[[x]]$xl/scale.factor, outlines[[x]]$yl/scale.factor, registration), border='#56565650', col=paste(as.character(registration$atlas$col[x]), '10', sep='') )})
        
index1<-which(!is.na(dataset1$acronym))
index2<-which(!is.na(dataset3$acronym))

section<-stereotactic.coordinates(outlines[[1]]$xl/scale.factor, outlines[[1]]$yl/scale.factor, registration)


axis(1, at=seq(-4,4,by=0.1), labels=FALSE, tck=-0.01, col.ticks='lightblue')
axis(1, at=seq(-4,4,by=0.5), labels=FALSE, tck=-0.015, col.ticks='coral')
axis(1, at=c(-4:4), labels=c(-4:4))

axis(2, at=seq(0,-6,by=-0.1), labels=FALSE, tck=-0.01, col.ticks='lightblue')
axis(2, at=seq(0,-6,by=-0.5), labels=FALSE, tck=-0.015, col.ticks='coral')
axis(2, at=c(0:-6), labels=c(0:-6), las=1)
box()


index1<-point.in.polygon(dataset1$ML, dataset1$DV, section$x, section$y)
index1<-which((!is.na(dataset1$acronym))&(index1))
index2<-point.in.polygon(dataset3$ML, dataset3$DV, section$x, section$y)
index2<-which((!is.na(dataset3$acronym))&(index2))
#index1<-1:nrow(dataset1)
#index2<-1:nrow(dataset3)

#symbols(dataset$ML[index1], dataset$DV[index1], circles=rep(0.05, length(index1)), inches=FALSE, bg='#fdb86395', asp=1,  add=T)

scale1<-data1[,which(names(data1)==ensemb)]
scale1 <-log2(scale1 +1)
scale1 <-(1-scale1/max(scale1) )

nlevels<-20

levels <- seq(0,1,length.out = nlevels)
col <- colorRampPalette(c("red","yellow"))(nlevels)  
colz <- col[cut(scale1,nlevels)]  


symbols(dataset1$ML[index1], dataset1$DV[index1], circles=rep(0.05, length(index1)), inches=FALSE, bg= colz, asp=1,  add=T)



scale2<-data2[,which(names(data2)==ensemb)]
scale2<-log2(scale2+1)
scale2<-(1-scale2/max(scale2) )

colz <- col[cut(scale2,nlevels)]  

symbols(dataset3$ML[index2], dataset3$DV[index2], circles=rep(0.05, length(index2)), inches=FALSE, bg= colz, asp=1,  add=T)

  

}