#!/usr/bin/Rscript

## Step 8: Project data in reduced dimensional space and generate UMAP plots and heatmaps

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Dimensionality reduction of MulTI-Tag data with SVD and UMAP
           
      Arguments:
      --k27=someValue    - H3K27me3 mapped fragment-to-peak file from Step 7
      --k4=someValue     - H3K4me2 mapped fragment-to-peak file from Step 7
      --k36=someValue    - H3K36me3 mapped fragment-to-peak file from Step 7
      --output=someValue - Output prefix

")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
invis <- gc(verbose=FALSE)

## Arg1 default
#if(is.null(args[1])){
if(is.null(argsL$input) | is.null(argsL$output)) {
  stop("Argument is missing!
      Dimensionality reduction of MulTI-Tag data with SVD and UMAP
           
      Arguments:
      --k27=someValue    - H3K27me3 mapped fragment-to-peak file from Step 7
      --k4=someValue     - H3K4me2 mapped fragment-to-peak file from Step 7
      --k36=someValue    - H3K36me3 mapped fragment-to-peak file from Step 7
      --output=someValue - Output prefix

")

  q(save="no")
}

library(ggplot2)

## Read in H3K27me3 data and cast into peak-by-cell matrix ##

temp<-read.table(argsL$k27)
temp2<-data.frame(paste(temp$V1,temp$V2,temp$V3, sep="_"), temp$V4, temp$V6)
colnames(temp2)<-c("feature", "cell", "count")
library(reshape)
temp3<-cast(temp2, formula= feature ~ cell)
temp3[is.na(temp3)]<-0
counts<-temp3[,c(2:ncol(temp3))]

## Scale fragment depth to median counts across cells ##

counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(median(apply(counts,2,sum))/apply(counts,2,sum))), nrow=nrow(counts)))
colnames(counts.termnorm)<-colnames(counts)

## Filter for top 40% of peaks by percentage of cells with fragments mapped ##

counts.termnorm.cellfilter<-counts.termnorm[apply(counts.termnorm, 1, function(x) length(x[x>0])) > quantile(apply(counts.termnorm, 1, function(x) length(x[x>0])), c(0.6)),]

## Perform TF-IDF normalization ##
                                                                                                                   
tfidf<-log(apply(counts.termnorm.cellfilter, 1, function(x) length(x)/length(x[x>0])))
counts.termnorm.cellfilter.tfidf<-counts.termnorm.cellfilter*tfidf
                 
## Log transform matrix ##
                 
counts.termnorm.cellfilter.tfidf.log<-log(counts.termnorm.cellfilter.tfidf+1)

## Perform Singular Value Decomposition and select components exceeding 0.2% of sum of diagonal terms ##

counts.termnorm.cellfilter.tfidf.log.svd<-svd(counts.termnorm.cellfilter.tfidf.log)
cutoff<-length(which(counts.termnorm.cellfilter.tfidf.log.svd$d/sum(counts.termnorm.cellfilter.tfidf.log.svd$d) > 0.002))
B<-counts.termnorm.cellfilter.tfidf.log.svd$u[,c(1:cutoff)]%*%diag(counts.termnorm.cellfilter.tfidf.log.svd$d[1:cutoff])
counts.termnorm.cellfilter.tfidf.log.svd.transform<-as.data.frame(t(B)%*%as.matrix(counts.termnorm.cellfilter.tfidf.log))
  
## Use UMAP to reduce dimensional space ##
 
library(umap)
counts.termnorm.cellfilter.tfidf.log.umap<-umap(t(counts.termnorm.cellfilter.tfidf.log.svd.transform))
frame<-as.data.frame(counts.termnorm.cellfilter.tfidf.log.umap$layout)
coldata<-unique(temp[,c(4:5)])
coldata<-coldata[order(coldata$V4),]
coldata2<-coldata[coldata$V4 %in% rownames(frame),]
frame$Source<-coldata2$V5
                 
## Plot data ##                 

k27.plot<-ggplot(frame, aes(x=V1, y=V2, color=Source)) + geom_point(cex=0.75) + theme_light() + xlab("UMAP 1") + ylab("UMAP2") + scale_color_manual(values=c("#F56647", "#1EA910")) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14))
#ggsave(paste(argsL$output, ".k27.UMAP.pdf", sep=""), k27.plot, width=8, height=7)

## Reassign counts for future merge with other target matrices ##

k27.counts<-counts


## Read in H3K4me2 data and cast into peak-by-cell matrix ##
                 
temp<-read.table(argsL$k4)
temp2<-data.frame(paste(temp$V1,temp$V2,temp$V3, sep="_"), temp$V4, temp$V6)
colnames(temp2)<-c("feature", "cell", "count")
library(reshape)
temp3<-cast(temp2, formula= feature ~ cell)
temp3[is.na(temp3)]<-0
counts<-temp3[,c(2:ncol(temp3))]
                 
 ## Perform same steps as above ##
                 
counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(median(apply(counts,2,sum))/apply(counts,2,sum))), nrow=nrow(counts)))
#counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(10000/apply(counts,2,sum))), nrow=nrow(counts)))
colnames(counts.termnorm)<-colnames(counts)
counts.termnorm.cellfilter<-counts.termnorm[apply(counts.termnorm, 1, function(x) length(x[x>0])) > quantile(apply(counts.termnorm, 1, function(x) length(x[x>0])), c(0.6)),]
tfidf<-log(apply(counts.termnorm.cellfilter, 1, function(x) length(x)/length(x[x>0])))
counts.termnorm.cellfilter.tfidf<-counts.termnorm.cellfilter*tfidf
counts.termnorm.cellfilter.tfidf.log<-log(counts.termnorm.cellfilter.tfidf+1)
counts.termnorm.cellfilter.tfidf.log.svd<-svd(counts.termnorm.cellfilter.tfidf.log)
cutoff<-length(which(counts.termnorm.cellfilter.tfidf.log.svd$d/sum(counts.termnorm.cellfilter.tfidf.log.svd$d) > 0.002))
B<-counts.termnorm.cellfilter.tfidf.log.svd$u[,c(1:cutoff)]%*%diag(counts.termnorm.cellfilter.tfidf.log.svd$d[1:cutoff])
counts.termnorm.cellfilter.tfidf.log.svd.transform<-as.data.frame(t(B)%*%as.matrix(counts.termnorm.cellfilter.tfidf.log))
library(umap)
counts.termnorm.cellfilter.tfidf.log.umap<-umap(t(counts.termnorm.cellfilter.tfidf.log.svd.transform))
frame<-as.data.frame(counts.termnorm.cellfilter.tfidf.log.umap$layout)
coldata<-unique(temp[,c(4:5)])
coldata<-coldata[order(coldata$V4),]
coldata2<-coldata[coldata$V4 %in% rownames(frame),]
frame$Source<-coldata2$V5
k4.plot<-ggplot(frame, aes(x=V1, y=V2, color=Source)) + geom_point(cex=0.75) + theme_light() + xlab("UMAP 1") + ylab("UMAP2") + scale_color_manual(values=c("#F56647", "#1EA910")) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14))
#ggsave(paste(argsL$output, ".k4.UMAP.pdf", sep=""), k4.plot, width=8, height=7)
k4.counts<-counts


## Read in H3K36me3 data and cast into peak-by-cell matrix ##
                 
temp<-read.table(argsL$k36)
temp2<-data.frame(paste(temp$V1,temp$V2,temp$V3, sep="_"), temp$V4, temp$V6)
colnames(temp2)<-c("feature", "cell", "count")
library(reshape)
temp3<-cast(temp2, formula= feature ~ cell)
temp3[is.na(temp3)]<-0
counts<-temp3[,c(2:ncol(temp3))]
                 
## Perform same steps as above ##   
                 
counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(median(apply(counts,2,sum))/apply(counts,2,sum))), nrow=nrow(counts)))
#counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(10000/apply(counts,2,sum))), nrow=nrow(counts)))
colnames(counts.termnorm)<-colnames(counts)
counts.termnorm.cellfilter<-counts.termnorm[apply(counts.termnorm, 1, function(x) length(x[x>0])) > quantile(apply(counts.termnorm, 1, function(x) length(x[x>0])), c(0.6)),]
tfidf<-log(apply(counts.termnorm.cellfilter, 1, function(x) length(x)/length(x[x>0])))
counts.termnorm.cellfilter.tfidf<-counts.termnorm.cellfilter*tfidf
counts.termnorm.cellfilter.tfidf.log<-log(counts.termnorm.cellfilter.tfidf+1)
counts.termnorm.cellfilter.tfidf.log.svd<-svd(counts.termnorm.cellfilter.tfidf.log)
cutoff<-length(which(counts.termnorm.cellfilter.tfidf.log.svd$d/sum(counts.termnorm.cellfilter.tfidf.log.svd$d) > 0.002))
B<-counts.termnorm.cellfilter.tfidf.log.svd$u[,c(1:cutoff)]%*%diag(counts.termnorm.cellfilter.tfidf.log.svd$d[1:cutoff])
counts.termnorm.cellfilter.tfidf.log.svd.transform<-as.data.frame(t(B)%*%as.matrix(counts.termnorm.cellfilter.tfidf.log))
library(umap)
counts.termnorm.cellfilter.tfidf.log.umap<-umap(t(counts.termnorm.cellfilter.tfidf.log.svd.transform))
frame<-as.data.frame(counts.termnorm.cellfilter.tfidf.log.umap$layout)
coldata<-unique(temp[,c(4:5)])
coldata<-coldata[order(coldata$V4),]
coldata2<-coldata[coldata$V4 %in% rownames(frame),]
frame$Source<-coldata2$V5
k36.plot<-ggplot(frame, aes(x=V1, y=V2, color=Source)) + geom_point(cex=0.75) + theme_light() + xlab("UMAP 1") + ylab("UMAP2") + scale_color_manual(values=c("#F56647", "#1EA910")) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14))
#ggsave(paste(argsL$output, ".k36.UMAP.pdf", sep=""), k36.plot, width=8, height=7)
k36.counts<-counts

## Merge all count matrices ##
                 
newcolnames<-substr(colnames(k27.counts), 1, 24)
colnames(k27.counts)<-newcolnames
colnames(k4.counts)<-newcolnames
colnames(k36.counts)<-newcolnames
counts<-rbind(k27.counts, k4.counts, k36.counts)
                 
## Perform same steps as above ##
                 
counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(median(apply(counts,2,sum))/apply(counts,2,sum))), nrow=nrow(counts)))
#counts.termnorm<-as.data.frame(as.matrix(t(t(counts)*(10000/apply(counts,2,sum))), nrow=nrow(counts)))
colnames(counts.termnorm)<-colnames(counts)
counts.termnorm.cellfilter<-counts.termnorm[apply(counts.termnorm, 1, function(x) length(x[x>0])) > quantile(apply(counts.termnorm, 1, function(x) length(x[x>0])), c(0.6)),]
tfidf<-log(apply(counts.termnorm.cellfilter, 1, function(x) length(x)/length(x[x>0])))
counts.termnorm.cellfilter.tfidf<-counts.termnorm.cellfilter*tfidf
counts.termnorm.cellfilter.tfidf.log<-log(counts.termnorm.cellfilter.tfidf+1)
counts.termnorm.cellfilter.tfidf.log.svd<-svd(counts.termnorm.cellfilter.tfidf.log)
cutoff<-length(which(counts.termnorm.cellfilter.tfidf.log.svd$d/sum(counts.termnorm.cellfilter.tfidf.log.svd$d) > 0.002))
B<-counts.termnorm.cellfilter.tfidf.log.svd$u[,c(1:cutoff)]%*%diag(counts.termnorm.cellfilter.tfidf.log.svd$d[1:cutoff])
counts.termnorm.cellfilter.tfidf.log.svd.transform<-as.data.frame(t(B)%*%as.matrix(counts.termnorm.cellfilter.tfidf.log))
library(umap)
counts.termnorm.cellfilter.tfidf.log.umap<-umap(t(counts.termnorm.cellfilter.tfidf.log.svd.transform))
frame<-as.data.frame(counts.termnorm.cellfilter.tfidf.log.umap$layout)
coldata<-unique(temp[,c(4:5)])
coldata<-coldata[order(coldata$V4),]
coldata$V4<-newcolnames
coldata2<-coldata[coldata$V4 %in% rownames(frame),]
frame$Source<-coldata2$V5
plot<-ggplot(frame, aes(x=V1, y=V2, color=Source)) + geom_point(cex=0.75) + theme_light() + xlab("UMAP 1") + ylab("UMAP2") + scale_color_manual(values=c("#F56647", "#1EA910")) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14))
#ggsave(paste(argsL$output, ".alltargets.UMAP.pdf", sep=""), plot, width=8, height=7)


## Select most variable peaks in first two components of SVD ##

comp1<-which(abs(counts.termnorm.cellfilter.tfidf.log.svd$u[,1]) > quantile(abs(counts.termnorm.cellfilter.tfidf.log.svd$u[,1]), 0.99))
comp2<-which(abs(counts.termnorm.cellfilter.tfidf.log.svd$u[,2]) > quantile(abs(counts.termnorm.cellfilter.tfidf.log.svd$u[,2]), 0.99))
comp<-unique(c(comp1,comp2))
regions<-data.frame(type=c(rep("K27me3", 4788), rep("K4me2", 15436), rep("K36me3", 21308)), color=c(rep("red", 4788), rep("purple", 15436), rep("#008080", 21308)))
RowCols<-regions$color[rownames(regions) %in% comp]

## Extract target-specific features from SVD and select most variable in first two components ## 

U<-as.data.frame(counts.termnorm.cellfilter.tfidf.log.svd$u)
rownames(U)<-rownames(counts.termnorm.cellfilter.tfidf.log)
U.k27<-U[as.numeric(rownames(U)) <= 4788,]
U.k4<-U[as.numeric(rownames(U)) > 4788 & as.numeric(rownames(U)) <= 20224,]
U.k36<-U[as.numeric(rownames(U)) > 20224,]
k27.comp1<-which(abs(U.k27[,1]) > quantile(abs(U.k27[,1]), 0.99))
k27.comp2<-which(abs(U.k27[,2]) > quantile(abs(U.k27[,2]), 0.99))
k4.comp1<-which(abs(U.k4[,1]) > quantile(abs(U.k4[,1]), 0.99))
k4.comp2<-which(abs(U.k4[,2]) > quantile(abs(U.k4[,2]), 0.99))
k36.comp1<-which(abs(U.k36[,1]) > quantile(abs(U.k36[,1]), 0.99))
k36.comp2<-which(abs(U.k36[,2]) > quantile(abs(U.k36[,2]), 0.99))
comp.mods<-unique(c(rownames(U.k27[k27.comp1,]),rownames(U.k27[k27.comp2,]),rownames(U.k4[k4.comp1,]),rownames(U.k4[k4.comp2,]),rownames(U.k36[k36.comp1,]),rownames(U.k36[k36.comp2,])))
RowCols.mods<-regions$color[rownames(regions) %in% comp.mods]

## Plot heatmaps ##

heatmap(as.matrix(counts.termnorm.cellfilter.tfidf.log[comp,]), col=my_palette, ColSideColors=collist, RowSideColors=RowCols)
heatmap(as.matrix(counts.termnorm.cellfilter.tfidf.log[comp.mods,]), col=my_palette, ColSideColors=collist, RowSideColors=RowCols.mods)
