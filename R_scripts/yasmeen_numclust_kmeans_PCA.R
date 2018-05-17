# This line is a comment! It starts with '#'
# infilename is the name of the input file. This should have only one header line
# and no lines above it. The first column should be called "Name"
infilename="E:/yasmeen/gp1_03222018.csv"
data=read.csv(infilename)
row.names(data)<-data$ï..Name 
data_mx=data.matrix(data)
data_mx <- data_mx[,-c(1)]
# Determine number of clusters
wss <- (nrow(data_mx)-1)*sum(apply(data_mx,2,var))
for (i in 2:25) wss[i] <- sum(kmeans(data_mx,centers=i)$withinss)
plot(1:25, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
# From scree plot elbow occurs at k = 8
# However, for simplicity we are just looking at k = 4 for now
# Apply k-means with k=3
k <- kmeans(data_mx, 4, nstart=25, iter.max=1000)
########## Run up to here for k-means clustering


# Scale
data2 <- data.frame(scale(data_mx))
# Verify variance is uniform
plot(sapply(data2, var)) 
pc <- prcomp(data2)
plot(pc)
summary(pc)
########## Run up to here to see variance covered by principal components..
########## ..looks as if the 5th PC is a natural elbow in the screen plot 
########## Also look at R console output (gives a more precise summary of cum variance)...
########## ..looks as if the first 5 PCs cover about 70% of the variance (70% or 90% are typical rules-of-thumb)


library(rgl)
comp <- data.frame(pc$x[,1:5])
# plot3d(comp$PC1, comp$PC2, comp$PC3, comp$PC4, comp$PC5)
plot3d(comp$PC1, comp$PC2, comp$PC3, col=k$clust, size=20)
########## Run up to here to see  principal components 1-3 of clusters
# plot3d(comp$PC1, comp$PC4, comp$PC5, col=k$clust, size=20)
########## Run up to here to see  principal components 1, 4, 5 of clusters
# Colors are in the following order, by cluster: black, red, green, blue, cyan, magenta... 
# For a complete list, visit: https://data.library.virginia.edu/setting-up-color-palettes-in-r/