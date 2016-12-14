#' Read a spatial transciptomics dataset.
#'
#' A TSV file that you want to read into working space and get collumn 1 and 2 as X and Y coordinates.
#' @filepath full filepath.
#' @examples
#' data.tmp<-read.ST(filepath='adjusted_stdata.tsv')
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

#' Symbolize EMSEMBLER IDs to REFseq symbols.
#'
#' bla bla bla.
#' @dataset data set read with read.ST.
#' @examples
#' data<-symbolize(data.tmp)
symbolize<-function(dataset){
new.names<-sapply(strsplit(names(dataset)[3:ncol(dataset)], '[.]'), '[', 1)

dataset<-dataset[,-(which(!(new.names%in%grcm38$ensgene))+2)]

new.names<-sapply(strsplit(names(dataset)[3:ncol(dataset)], '[.]'), '[', 1)

symbol<-grcm38$symbol[-which(duplicated(grcm38$ensgene))]
ensgene<-grcm38$ensgene[-which(duplicated(grcm38$ensgene))]

#replace with symbol
new.names<-symbol[which(ensgene%in%new.names)]

dataset<-dataset[, -which(duplicated(new.names))]
names(dataset)[3:ncol(dataset)]<-new.names[-which(duplicated(new.names))]

return(dataset)
}


#EXAMPLES

data.tmp<-read.ST('~/GitHub/St-dat/Sequencing_reads/E2_ID3_S3_stdata_adj.csv')
data<-symbolize(data.tmp)