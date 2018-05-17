df2 <- data.frame(supp=rep(rep(c("tcga", "msk"),each=2),2),
			traits=c(rep(c("Male", "Female"), 2),rep(c("Hepatitis B", "Hepatitis C"),2)),
                  len=c(129, 67, 103, 33, 44, 35, 22, 50))
p <- ggplot(data=df2, aes(supp, y=len, fill=traits)) +
  geom_bar(stat="identity", position="stack", width=0.05)

p + coord_flip()