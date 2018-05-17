## test <- data.frame(institute=c("MSK","TCGA"),val1=c(22,16),val2=c(18,37),val3=c(60,55),
## val4=c(34,23),val5=c(66,77))

test <- data.frame(variable=c("Viral status","Gender"),hepb_tcga=c(22,0),hepc_tcga=c(18,0),novir_tcga=c(60,0),female_tcga=c(0,34),male_tcga=(0,66),hepb_msk=c(16,0),hepc_msk=c(37,0),novir_msk=c(53,0),female_msk=c(0,23),male_msk=(0,77))

test  <- data.frame(variable=c("Viral status","Gender"), 
                    hepb_tcga=c(22,0),     
                    hepc_tcga=c(18,0), 
			  novir_tcga=c(60,0),   
			  female_tcga=c(0,34),   
			  male_tcga=c(0,66),   
			  hepb_msk=c(16,0),   
			  hepc_msk=c(37,0),   
			  novir_msk=c(53,0),   
			  female_msk=c(0,23),   
                    male_msk=c(0,77))

library(reshape2) 
library(ggplot2)
melted <- melt(test, "variable")
melted$cat <- ''
melted[melted$variable == 'hepb_tcga',]$cat <- "TCGA"
melted[melted$variable == 'hepc_tcga',]$cat <- "TCGA"
melted[melted$variable == 'novir_tcga',]$cat <- "TCGA"
melted[melted$variable == 'female_tcga',]$cat <- "TCGA"
melted[melted$variable == 'male_tcga',]$cat <- "TCGA"
melted[melted$variable == 'hepb_msk',]$cat <- "MSK"
melted[melted$variable == 'hepc_msk',]$cat <- "MSK"
melted[melted$variable == 'novir_msk',]$cat <- "MSK"
melted[melted$variable == 'female_msk',]$cat <- "MSK"
melted[melted$variable == 'male_msk',]$cat <- "MSK"

ggplot(melted, aes(x = cat, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ institute)