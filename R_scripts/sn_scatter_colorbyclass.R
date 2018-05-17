infilename = "E:/Subiksha/Gene frequency, mutation count etc. - Sheet16 (2).csv"
hm=read.csv(infilename)
hm=data.frame(idsupp=hm$Sample_ID ,mutationcount=hm$Mutation_Count,cna=hm$CNA,group =hm$Group,dat=hm$Z)


c1 <- plot(hm[,c(3,2)], col= c("red","blue") [hm$dat])
points(hm[which(hm$group=='Selected'),c(3,2)], col= c("red","blue") [hm$dat])
