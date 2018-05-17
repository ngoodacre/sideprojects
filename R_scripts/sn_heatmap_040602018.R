# This line is a comment! It starts with '#'
# infilename is the name of the input file. This should have only one header line
# and no lines above it. The first column should be called "Name"
infilename="E:/S/test2.csv"
hm=read.csv(infilename)
row.names(hm)<-hm$ï..Name 
hm_mx=data.matrix(hm)
hm_mx <- hm_mx[,-c(1)]
# the scale() function normalizes linearly from -1 to 1. There are other ways to normalize, e.g. sigmoid
# We can look into these as well
hm_mx_scale=scale(hm_mx)
# play around with the n value.. a higher number gives a smoother (but noisier) color distribution
# while a smaller number gives a simpler, but grainier look
my_palette <- colorRampPalette(c("red","yellow","green"))(n = 10)
# The Colv=NA parameter removes the column dendrogram
# Leave Colv paramter out to include dendrogram for columns by default
# the margins parameter increases the margin width in the image. 
# heatmap(hm_mx_scale,distfun=corr.dist, hclustfun=avg, col=my_palette,margins=c(10,10))
hc <- hclust(as.dist(1-cor(t(hm_mx_scale))))