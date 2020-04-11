\name{BioMGC}
\alias{BioMGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extraction of Biomarker Co-clusters using Individual Control Chart
}
\description{
Function for Calculating Co-cluster Means, Ranking Them and Their Corresponding Gene and Compound Clusters,Extraction of Biomarker Co-clusters using using Individual Control Chart
}
\usage{
BioMGC(data,rk,ck)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rk}{
row clusters estimated observing dendrogram
}
\item{ck}{
column clusters estimated observing dendrogram
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \link{list} of object that containing the following:
\itemize{
\item \code{NGC_cls_MeanMat}: Ranked Co-Cluster Mean and Their Corresponding Row and Column Cluster. The First number of column "NGC_Cocls" is row cluster and the second number is the column cluster.
\item \code{CoClsDtMat}: Re-organized Data matrix based on "rHC_Function"
\item \code{Co_Gcls}: Gene clusters.
\item \code{Co_Ccls}: Column clusters.
\item \code{colorsG}: Color of gene (row) cluster for generating co_cluster graph.
\item \code{colorsC}: Color of Compound (column) cluster for generating co_cluster graph.
\item \code{CentralLine}: Central line of QCC for identification of Biomarker Co-Cluster.
\item \code{UpContLimit}: Uper Control Limit (UCL) of QCC for identification of Biomarker Co-Cluster.
\item \code{LowrContLimit}: Lower Control Limit (LCL) of QCC for identification of Biomarker Co-Cluster.
\item \code{color}: Color for QCC.
\item \code{pchmark}: PCmark for QCC.
}
}

\author{
Md Bahadur Badsha
}

\examples{
###### Real Data 
dataExpr <- as.matrix(read.csv(file="GluDtTime2.csv", row.names=1))
#dataExpr <- as.matrix(read.csv(file="PPARdatAllTime.csv", row.names=1))
GCmat <- dataExpr
##------------------------------------------------------------------------

###### Simulated Data 
GCmat <- as.matrix(read.csv(file="SimulationData.csv", row.names=1))
##------------------------------------------------------------------------

###### Data Analysis Code using "rHC_Function" 
source("rHC_Function.R")
BioGnCm <- BioMGC(GCmat, 3, 3)  #names(BioGnCm)

## Ranked Co-Cluster Mean and Their Corresponding Row and Column Cluster 
## The First number of column "NGC_Cocls" is row cluster and the second number 
## is the column cluster  
GC_cls_MeanMat <- BioGnCm$NGC_cls_MeanMat

## Re-organized Data matrix based on "rHC_Function"
CoClsDtMat <- BioGnCm$CoClsDtMat

## Gene clusters 
Co_Gcls <- BioGnCm$Co_Gcls

## Column clusters 
Co_Ccls <- BioGnCm$Co_Ccls

## Color of gene (row) cluster for generating co_cluster graph 
colorsG <- BioGnCm$colorsG

## Color of Compound (column) cluster for generating co_cluster graph 
colorsC <- BioGnCm$colorsC

## Central line of QCC for identification of Biomarker Co-Cluster 
CntrLine_QC <- BioGnCm$CentralLine

## Uper Control Limit (UCL) of QCC for identification of Biomarker Co-Cluster 
UCL_QC <- BioGnCm$UpContLimit

## Lower Control Limit (LCL) of QCC for identification of Biomarker Co-Cluster
LCL_QC <- BioGnCm$LowrContLimit

## Color QCC
ColorQC <- BioGnCm$color

## PCmark for QCC
PcmQC <- BioGnCm$pchmark
##--------------------------------------------------------------------------------------------------------

## Gene (Row) and Compound (Column) Co-Cluster Graph 
library(fields)
par(mar=c(4,8,1,5.5))
image(CoClsDtMat, col=colorRampPalette(c("white","gray","black"), space = "rgb")(300),axes=FALSE)
mtext(text=rownames(CoClsDtMat), side=1, line=0.3, at=seq(0,1,1/(nrow(CoClsDtMat)-1)), las=2, cex=0.6, col=colorsG)
mtext(text=colnames(CoClsDtMat), side=2, line=0.3, at=seq(0,1,1/(ncol(CoClsDtMat)-1)), las=2, cex=0.6, col=colorsC)
image.plot(CoClsDtMat, col=colorRampPalette(c("white","gray","black"), space = "rgb")(300),legend.only=TRUE, horizontal=F)
##--------------------------------------------------------------------------------------------------------

## Graph of QCC for identification of biomarker co-cluster
plot(GC_cls_MeanMat[,2],xlab="Combination of Row and Column Cluster",ylab="Co-cluster Average",xaxt='n',
	pch=PcmQC,cex=1.1,col=ColorQC, ylim = c(min(LCL_QC, min(GC_cls_MeanMat[,2]-5)), max(GC_cls_MeanMat[,2]+2)))
abline(h=c(UCL_QC, CntrLine_QC, LCL_QC), lty=3, lwd=3, col="green")
mtext(text = GC_cls_MeanMat[,1], side=1, line=0.3, at=seq(1:length(GC_cls_MeanMat[,1])), las=2, font=2,
	cex=1, col="black")
text(c(length(GC_cls_MeanMat[,1])-0.5,length(GC_cls_MeanMat[,1])-0.5,length(GC_cls_MeanMat[,1])-0.5),
	c(UCL_QC, CntrLine_QC, LCL_QC), c("UCL","CL","LCL"),font=4,cex=1.5,col="blue")
##---------------------------------------------------------------------------------------------------------
}

