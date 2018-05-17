# This line is a comment! It starts with '#'
# infilename is the name of the input file. This should have only one header line
# and no lines above it. The first column should be called "Name"
infilename="E:/yasmeen/gp1_pc1.csv"
hm=read.csv(infilename)
row.names(hm)<-hm$Name
hm_mx=data.matrix(hm)
x=t(hm_mx)
hm_mx <- hm_mx[,-c(1)]
x <- x[-c(1),]
# play around with the n value.. a higher number gives a smoother (but noisier) color distribution
# while a smaller number gives a simpler, but grainier look
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)
library(RColorBrewer)
hmcols<- colorRampPalette(brewer.pal(9,"GnBu"))(100)
# cols = as.numeric(as.factor(colnames(x)))
library(gplots)
# This is where you label the columns
v=c(rep("Primary Site",10),rep("Metastasis",21))
# This is where you color according to label of column
# cols<- sample(colors(), 31, replace = F)[as.numeric(as.factor((v)))]
cols<- c("red","darkorange1")[as.numeric(as.factor((v)))]
# cols=c(rep("darkorange1",10),rep("red",21))
## The variables to explore here are: cor(method=?), hclustfun=function(x)hclust(c,method=?), scale=?
## cor(method=?) is the correlation method. Standard is "pearson", other options include: "kendall", "spearman"
## hclust(c,method=?) is the clustering method. Standard is "ward.D", other options: "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
## scale=? is the normalization method. Standard is "row"
## The Rowv=NA parameter removes the row dendrogram
heatmap(x, distfun=function(x) as.dist(1-cor(t(x),method="pearson")), 
hclustfun=function(x)hclust(x,method="ward.D"),trace="none", scale = "row", 
col=my_palette, labRow = rownames(x), Rowv=NA, margins = c(10,10), density.info="none",cexRow=1.3,cexCol=1.0)

