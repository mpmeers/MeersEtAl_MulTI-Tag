#!/usr/bin/Rscript

## Step 5: Filter cells based on unique fragment counts across all targets

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 10) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Filter MulTI-Tag bed files based on unique reads
      
      Arguments:
      --k27=someValue    - H3K27me3 CellRanger-formatted bed file
      --k4=someValue     - H3K4me2 CellRanger-formatted bed file
      --k36=someValue    - H3K36me3 CellRanger-formatted bed file
      --k27_H1=somevalue - H3K27me3 unique fragments in H1 cells text file
      --k4_H1=somevalue  - H3K4me2 unique fragments in H1 cells text file
      --k36_H1=somevalue - H3K36me3 unique fragments in H1 cells text file
      --k27_K5=somevalue - H3K27me3 unique fragments in K562 cells text file
      --k4_K5=somevalue  - H3K4me2 unique fragments in K562 cells text file
      --k36_K5=somevalue - H3K36me3 unique fragments in K562 cells text file
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
      Filter MulTI-Tag bed files based on unique reads
      
      Arguments:
      --k27=someValue    - H3K27me3 CellRanger-formatted bed file
      --k4=someValue     - H3K4me2 CellRanger-formatted bed file
      --k36=someValue    - H3K36me3 CellRanger-formatted bed file
      --k27_H1=somevalue - H3K27me3 unique fragments in H1 cells text file
      --k4_H1=somevalue  - H3K4me2 unique fragments in H1 cells text file
      --k36_H1=somevalue - H3K36me3 unique fragments in H1 cells text file
      --k27_K5=somevalue - H3K27me3 unique fragments in K562 cells text file
      --k4_K5=somevalue  - H3K4me2 unique fragments in K562 cells text file
      --k36_K5=somevalue - H3K36me3 unique fragments in K562 cells text file
      --output=someValue - Output prefix

")

  q(save="no")
}

k27<-rbind(read.table(argsL$k27_H1), read.table(argsL$k27_K5))
k4<-rbind(read.table(argsL$k27_H1), read.table(argsL$k27_K5))
k36<-rbind(read.table(argsL$k27_H1), read.table(argsL$k27_K5))
colnames(k27)<-c("Fragments", "Barcode")
colnames(k4)<-c("Fragments", "Barcode")
colnames(k36)<-c("Fragments", "Barcode")
k27<-k27[order(k27$Fragments, decreasing=TRUE),]
k4<-k4[order(k4$Fragments, decreasing=TRUE),]
k36<-k36[order(k36$Fragments, decreasing=TRUE),]
k27$Rank<-c(1:nrow(k27))
k4$Rank<-c(1:nrow(k4))
k36$Rank<-c(1:nrow(k36))
k27$Target<-"H3K27me3"
k4$Target<-"H3K4me2"
k36$Target<-"H3K36me3"
#all<-rbind(k27,k4,k36)
#library(ggplot2)
#ggplot(all, aes(x=Rank, y=Fragments, color=Target)) + geom_line(lwd=1.5) + theme_light() + scale_x_log10() + scale_y_log10()

k27.bed<-read.table(argsL$k27)
k4.bed<-read.table(argsL$k4)
k36.bed<-read.table(argsL$k36)
k27.bed$V6<-substr(k27.bed$V4, 1, 24)
k4.bed$V6<-substr(k4.bed$V4, 1, 24)
k36.bed$V6<-substr(k36.bed$V4, 1, 24)
k27$Barcode.cell<-substr(k27$Barcode, 1, 24)
k4$Barcode.cell<-substr(k4$Barcode, 1, 24)
k36$Barcode.cell<-substr(k36$Barcode, 1, 24)
k27.bed2<-k27.bed[k27.bed$V6 %in% k27$Barcode.cell[k27$Fragments > 500] & k27.bed$V6 %in% k4$Barcode.cell[k4$Fragments > 200] & k27.bed$V6 %in% k36$Barcode.cell[k36$Fragments > 200],c(1:5)]
k4.bed2<-k4.bed[k4.bed$V6 %in% k27$Barcode.cell[k27$Fragments > 500] & k4.bed$V6 %in% k4$Barcode.cell[k4$Fragments > 200] & k4.bed$V6 %in% k36$Barcode.cell[k36$Fragments > 200],c(1:5)]
k36.bed2<-k36.bed[k36.bed$V6 %in% k27$Barcode.cell[k27$Fragments > 500] & k36.bed$V6 %in% k4$Barcode.cell[k4$Fragments > 200] & k36.bed$V6 %in% k36$Barcode.cell[k36$Fragments > 200],c(1:5)]
k27.k5<-read.table(argsL$k27_K5)
k27.h1<-read.table(argsL$k27_H1)
k4.k5<-read.table(argsL$k4_K5)
k4.h1<-read.table(argsL$k4_H1)
k36.k5<-read.table(argsL$k36_K5)
k36.h1<-read.table(argsL$k36_H1)
k27.bed2$V6<-"none"
k27.bed2$V6[k27.bed2$V4 %in% k27.h1$V2]<-"H1"
k27.bed2$V6[k27.bed2$V4 %in% k27.k5$V2]<-"K562"
k4.bed2$V6<-"none"
k4.bed2$V6[k4.bed2$V4 %in% k4.h1$V2]<-"H1"
k4.bed2$V6[k4.bed2$V4 %in% k4.k5$V2]<-"K562"
k36.bed2$V6<-"none"
k36.bed2$V6[k36.bed2$V4 %in% k36.h1$V2]<-"H1"
k36.bed2$V6[k36.bed2$V4 %in% k36.k5$V2]<-"K562"
write.table(k27.bed2, paste(argsL$output, ".k27.bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(k4.bed2, paste(argsL$output, ".k4.bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(k36.bed2, paste(argsL$output, ".k36.bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
