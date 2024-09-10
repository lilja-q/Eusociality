###R code used for analysis of diversification with MEDUSA in Blaimer et al.XXXX ###

library(ape)
library(geiger)

HymTopCf <- read.tree(file="TopCfamily.tre")
HymTopAf <- read.tree(file="TopAfamily.tre")

richness <- read.csv(file = "richness.csv", header = T) ## header needs to be specific format, check example

TopCres <- medusa(HymTopCf, richness = richness, criterion = "aicc", cut= "node", ncores=1)
TopAres <- medusa(HymTopAf, richness = richness, criterion = "aicc", cut= "node", ncores=1)


write.csv(TopCres[["zSummary"]], file="TopCres.csv")
write.csv(TopAres[["zSummary"]], file="TopAres.csv")
plot(TopAres, cex=0.5, time=TRUE)
plot(TopCres, cex=0.5, time=TRUE)