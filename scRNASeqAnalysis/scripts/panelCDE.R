library(entropy)
library(Rtsne)



samples<-c("01","02")

#read CellRanger 10x output
G<-read.table("./01/features.tsv",header=FALSE,sep="\t",stringsAsFactors=FALSE)
ng<-dim(G)[1]

nb2<-rep(0,2) #numbers of bar codes per sample

CC<-list()
for (i in 1:2) {
   sname<-samples[i]
   barcodes<-read.table(paste("./",sname,"/barcodes.tsv",sep=""),header=FALSE)
   M<-read.table(paste("./",sname,"/matrix.mtx",sep=""),header=FALSE,skip=3,sep=" ")
   nb<-dim(barcodes)[1]
   C<-matrix(0,nrow=ng,ncol=nb)
   C[cbind(M[,1],M[,2])]<-M[,3]
   rownames(C)<-G[,2] #gene symbols
   expressed<-rowSums(C)>0
   CC[[i]]<-C[expressed,]
   nb2[i]<-nb
}



#create a single matrix of counts from both samples
x1<-CC[[1]]
x2<-CC[[2]]

cnb<-cumsum(nb2)
genes<-union(rownames(x1),rownames(x2))
ng<-length(genes)
P<-matrix(0,nrow=ng,ncol=cnb[2])
rownames(P)<-genes

#put x1 into P
indx<-match(rownames(x1),genes)
P[indx,1:nb2[1]]<-x1
#put x2 into P
indx<-match(rownames(x2),genes)
P[indx,(cnb[1]+1):cnb[2]]<-x2
#P contains read counts; genes in rows, cells in columns

#filtering by mean expression count
detected<-rowMeans(P)>.1
P<-P[detected,] 
csP<-colSums(P)
#normalization; cells in rows, genes in columns:
P<-t(P)/csP
#P[i,j] is the probability that a randomly chosen read from cell i will be from gene j 




#tSNE work
N<-9000 #use the same number of cells from both samples, for balance
P<-P[c(sample(nb2[1],N),sample(nb2[2],N)+nb2[1]),]
n<-2*N

#D is the Jensen-Shannon distance matrix
D<-matrix(0,ncol=n,nrow=n)
E<-apply(P,1,entropy)
for (i in 1:(n-1)) {
   for (j in (i+1):n) {
      D[i,j]<-sqrt(entropy((P[i,]+P[j,])/2)-E[i]/2-E[j]/2)
   }
}
D<-D+t(D)

#tsne plots
set.seed(2357)
tsne<-Rtsne(D,is_distance=TRUE,perplexity=50,check_duplicates=FALSE)

pdf("panelC_tSNE.pdf")
plot(tsne$Y,pch="",main="sample 1, all cells",xlab="tSNE 1",ylab="tSNE 2")
points(tsne$Y[1:cnb[1],],col="blue",pch=16,cex=.3)

plot(tsne$Y,pch="",main="sample 2, all cells",xlab="tSNE 1",ylab="tSNE 2")
points(tsne$Y[(cnb[1]+1):cnb[2],],col="red",pch=16,cex=.3)
dev.off()




#tSNE plots based on a subset of genes, a "geneset" 
#distance must be calculated with account of genes NOT in the geneset as a single collection bin "binN"
comparisons<-c("AMLvAN","MFvAN")
for (comparison in comparisons) {
   #read a gene set
   GS<-read.table(paste0("./",comparison,"_geneset.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="\"")
   geneset<-GS[,1]
   ng<-length(geneset)

   in.set<-colnames(P) %in% geneset
   binN<-rowSums(P[,!in.set])
   xP<-cbind(P[,in.set],binN)
   #this is the probability matrix with the mandatory collection bin in the last row, which is the probability that a randomly chosen read is NOT from any gene in the geneset
   #rowSums(xP) == 1 is still true


   #Jensen-Shannon distance matrix:
   D<-matrix(0,ncol=n,nrow=n)
   E<-apply(xP,1,entropy)
   for (i in 1:(n-1)) {
      for (j in (i+1):n) {
         D[i,j]<-sqrt(entropy((xP[i,]+xP[j,])/2)-E[i]/2-E[j]/2)
      }
   }
   D<-D+t(D)

   set.seed(2357)
   tsne<-Rtsne(D,perplexity=50,is_distance=TRUE)
   
   pdf(paste0(comparison,"_tSNE.pdf"))
   plot(tsne$Y,pch=16,cex=.3,type="n",xlab="tSNE 1",ylab="tSNE 2")
   points(tsne$Y[1:cnb[1],],pch=16,cex=.3,col="blue")
   points(tsne$Y[(cnb[1]+1):cnb[2],],pch=16,cex=.3,col="red")
   dev.off()
}

