test1 <- read.csv("E:/Subiksha/MSK_TCGA/MSK_TCGA_tumordata_04252018.csv")
test=data.frame(supp=test1$ï..supp ,v=test1$variable,vtype=test1$variable_type,dat=test1$data)

#test  <- data.frame(person=c("A", "B", "C", "D", "E"), 
#                    value1=c(100,150,120,80,150),
#			  value2=c(25,30,45,30,30), 
#                    value3=c(100,120,150,150,200)) 

library(reshape2) # for melt

#melted <- melt(test, "supp")
#melted$cat <- ''
#melted[melted$variable == 'value1',]$cat <- "TCGA"
#melted[melted$variable != 'value1',]$cat <- "MSK"

p1 <- ggplot(test, aes(x = supp, y = dat, fill = vtype)) +  geom_bar(stat = 'identity', position = 'stack') + facet_grid(v ~ .)
p1 <- p1 + coord_flip()
p1

