#########################################################################
#######################	Multiple plot function	#########################
#########################################################################
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#########################################################################
####################	Execute code	###############################
#########################################################################

df1 <- data.frame(supp=rep(c("tcga", "msk"),each=2),
			gender=rep(c("Male", "Female"), 2),
                  len=c(129, 67, 103, 33))


df2 <- data.frame(supp=rep(c("tcga", "msk"),each=2),
			virus=rep(c("Hepatitis B", "Hepatitis C"),2),
			len=c(44, 35, 22, 50))

## df3=rbind(df1,df2)

par(mfrow=c(2,1))
p1 <- ggplot(data=df1, aes(supp, y=len, fill=gender)) +
  geom_bar(stat="identity", position="stack", width=0.05,xlab("supplier"),ylab("patients"))
p1 <- p1 + labs(x="supplier",y="patients")
p1 <- p1 + coord_flip()
p2 <- ggplot(data=df2, aes(supp, y=len, fill=virus)) +
  geom_bar(stat="identity", position="stack", width=0.05)
p2 <- p2 + labs(x="supplier",y="patients")
p2 <- p2 + coord_flip()
multiplot(p1, p2,rows=2)