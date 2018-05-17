# This line is a comment! It starts with '#'
# infilename is the name of the input file. This should have only one header line
# and no lines above it. The first column should be called "Name"
infilename="E:/Subiksha/test2.csv"
hm=read.csv(infilename)
row.names(hm)<-hm$Cancer 
hm_mx=data.matrix(hm)
hm_mx <- hm_mx[,-c(1)]
x=t(x)
rownames(x)=colnames(hm_mx)
colnames(x)=rownames(hm_mx)
# play around with the n value.. a higher number gives a smoother (but noisier) color distribution
# while a smaller number gives a simpler, but grainier look
my_palette <- colorRampPalette(c("blue","grey","red"))(n = 25)
library(RColorBrewer)
hmcols<- colorRampPalette(brewer.pal(9,"GnBu"))(100)
cols = as.numeric(as.factor(colnames(x)))
# The Colv=NA parameter removes the column dendrogram
heatmap.2(x, distfun=function(x) as.dist(1-cor(t(x),method="pearson")), 
hclustfun=function(x)hclust(x,method="ward.D"),trace="none", scale = "row", 
col=my_palette, labCol=colnames(x), labRow = rownames(x), margins = c(10,10), density.info="none", 
main = "scaled gene and correlation distance",cexRow=0.5,cexCol=0.75)
 