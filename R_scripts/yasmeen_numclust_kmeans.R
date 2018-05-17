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
for (i in 2:25) wss[i] <- sum(kmeans(data_mx,
                                     centers=i)$withinss)
plot(1:25, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
####### Run up to here first to see # clusters plot


# From scree plot elbow occurs at k = 10
# However, for simplicity we are just looking at k = 3 for now
# Apply k-means with k=3
k <- kmeans(data_mx, 3, nstart=25, iter.max=1000)
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
par(mfrow=c(3,3))
# i and j are columns. Vary (i in 1:25) and (j in 1:25)... for instance try (i in 1:3), (j in 1:3)
# then (i in 1:3), (j in 4:6), then (i in 4:6), (j in 1:3), etc. Save your output. 
# ultimately you want to look at ALL possible i,j combinations (there are 225)
for (i in 1:25){
	xdata=data_mx[,i]
	for (j in 1:25){
		ydata=data_mx[,j]
		plot(xdata,ydata,xlab=colnames(data_mx)[i],ylab=colnames(data_mx)[j],col=k$clust,pch=16)
		   		}
			}
####### Run up to here to see clustering across pairs of columns