df2 <- data.frame(supp=rep(c("tcga", "msk"),each=2),
			gender=rep(c("Male", "Female"), 2),
                  len=c(67, 129, 33, 42))
p <- ggplot(data=df2, aes(supp, y=len, fill=gender)) +
  geom_bar(stat="identity", position="stack", width=0.05)

p + coord_flip()