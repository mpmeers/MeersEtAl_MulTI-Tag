#!/usr/bin/Rscript

## Step 8: Project data in reduced dimensional space

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
temp<-read.table(argsL$k27)
temp2<-data.frame(paste(temp$V1,temp$V2,temp$V3, sep="_"), temp$V4, temp$V6)
colnames(temp2)<-c("feature", "cell", "count")
library(reshape)
temp3<-cast(temp2, formula= feature ~ cell)
temp3[is.na(temp3)]<-0
counts<-temp3[,c(2:ncol(temp3))]
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
k27.plot<-ggplot(frame, aes(x=V1, y=V2, color=Source)) + geom_point(cex=0.75) + theme_light() + xlab("UMAP 1") + ylab("UMAP2") + scale_color_manual(values=c("#F56647", "#1EA910")) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14))
#ggsave(paste(argsL$output, ".k27.UMAP.pdf", sep=""), k27.plot, width=8, height=7)
k27.counts<-counts

temp<-read.table(argsL$k4)
temp2<-data.frame(paste(temp$V1,temp$V2,temp$V3, sep="_"), temp$V4, temp$V6)
colnames(temp2)<-c("feature", "cell", "count")
library(reshape)
temp3<-cast(temp2, formula= feature ~ cell)
temp3[is.na(temp3)]<-0
counts<-temp3[,c(2:ncol(temp3))]
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

temp<-read.table(argsL$k36)
temp2<-data.frame(paste(temp$V1,temp$V2,temp$V3, sep="_"), temp$V4, temp$V6)
colnames(temp2)<-c("feature", "cell", "count")
library(reshape)
temp3<-cast(temp2, formula= feature ~ cell)
temp3[is.na(temp3)]<-0
counts<-temp3[,c(2:ncol(temp3))]
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

newcolnames<-substr(colnames(k27.counts), 1, 24)
colnames(k27.counts)<-newcolnames
colnames(k4.counts)<-newcolnames
colnames(k36.counts)<-newcolnames
counts<-rbind(k27.counts, k4.counts, k36.counts)
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
#ggsave("plots/MPM2244-2246.sci.MulTI.K27gt500-K4m1gt200-k36gt200uniquereads.bedgraph.stringent.singlecellmap.bed_UMAP.pdf", plot, width=8, height=7)
