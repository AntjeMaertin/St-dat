
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



Min <- pmin(dataset2$x, dataset2$y) 
Max <- pmax(dataset2$x, dataset2$y) 


corner <- Min != Max 

quartz()
plot(dataset2$x, dataset2$y, bg=dataset2$color, pch=21, ylim=c(regi$transformationgrid$height,0), xlim=c(0,regi$transformationgrid$width), asp=1)
library(geometry)
MBR <- function(points) {
    tryCatch({
        a2 <- geometry::convhulln(points, options = 'FA')

        e <- points[a2$hull[,2],] - points[a2$hull[,1],]            # Edge directions
        norms <- apply(e, 1, function(x) sqrt(x %*% x)) # Edge lengths

        v <- diag(1/norms) %*% as.matrix(e)                        # Unit edge directions


        w <- cbind(-v[,2], v[,1])                       # Normal directions to the edges

        # Find the MBR
        vertices <- as.matrix((points) [a2$hull, 1:2])    # Convex hull vertices
        minmax <- function(x) c(min(x), max(x))         # Computes min and max
        x <- apply(vertices %*% t(v), 2, minmax)        # Extremes along edges
        y <- apply(vertices %*% t(w), 2, minmax)        # Extremes normal to edges
        areas <- (y[1,]-y[2,])*(x[1,]-x[2,])            # Areas
        k <- which.min(areas)                           # Index of the best edge (smallest area)

        # Form a rectangle from the extremes of the best edge
        as.data.frame(cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,]))
    }, error = function(e) {
        assign('points', points, .GlobalEnv)
        stop(e)  
    })
}



array.corners<-MBR(cbind(dataset2$x, dataset2$y))

xcoord<-seq(array.corners[1,1], array.corners[4,1], length.out=35)
ycoord<-seq(array.corners[1,2], array.corners[2,2], length.out=33)

ycoord<-rep(ycoord, 35)
xcoord <-rep(xcoord, 33)

indexX<-rep(1:35, 33)

indexY<-rep(1:33, 35)

points(xcoord, ycoord, col='red')
text(xcoord, ycoord, labels=paste(indexX, indexY, sep='x'), cex=0.5)

distance<-sqrt( (array.corners[1,1]-dataset2$x)^2 + (array.corners[1,2]-dataset2$y)^2 )
corner <-which.min(distance)

points(dataset2$x[corner], dataset2$y[corner], pch=16, col='pink')

quartz()
plot(as.numeric(as.vector(dat1[,1])), as.numeric(as.vector(dat1[,2])), cex=1, ylim=c(32, 2), xlim=c(34,2))
points(as.numeric(as.vector(dat1.adj[,1])), as.numeric(as.vector(dat1.adj[,2])), cex=1, pch=21, bg='red')


names(dat1)[3:ncol(dat1)]%in%

