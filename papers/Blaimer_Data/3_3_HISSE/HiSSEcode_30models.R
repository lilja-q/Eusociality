#####Assessing trait-dependent diversification within the HiSSE framework####
# This code fits and evaluates various models Bisse and HiSSE models, including several null models
# Sources are:
# 1) the HiSSE vignette by Jeremy Beaulieu included in the HiSSE package (v1.9.18)
# 2) the code provided by Beaulieu and O'Meara (2016) and 
# 3) code published by Harrington & Reeder (2016)

#Load relevant packages
library(ape)
library(geiger)
library(diversitree)
library(hisse)
library(parallel)

#download HiSSE.null4.9rate.R from Dryad repository https://datadryad.org/stash/dataset/doi:10.5061/dryad.c494p
load(HiSSE.null4.9rate.R)

TopA <- read.tree(file="HYMtopA_bamm_hisse.phy")
TopC <- read.tree(file="HYMtopC_bamm_hisse.phy")
Sting <- read.csv("Sting.csv", header=FALSE)
Waist <- read.csv("ww.csv", header=FALSE)
Para <- read.csv("Parasitoid.csv", header=FALSE)
entomo <- read.csv("Entomophagy.csv", header=FALSE)
phyto <- read.csv("Phytophagy.csv", header=FALSE)
secphyto <- read.csv("SecPhytophagy.csv", header=FALSE)

sampfrac_s <-c(0.007354953, 0.001846296) # sampling frequency for stinger
sampfrac_w <-c(0.003679269, 0.005082557) # sampling frequency for waist
sampfrac_para <-c(0.002610357, 0.006081215) # sampling frequency for parasitoidism
sampfrac_entomo <-c(0.003443052, 0.005378214) # sampling frequency for entomophagy
sampfrac_phyto <-c(0.006019132, 0.003784913) # sampling frequency for phytophagy
sampfrac_secphyto <-c(0.005801694, 0.003825655) # sampling frequency for secondary phytophagy


##Set up all relevant transition matrices##

trans.rates.hisse <- TransMatMaker.old(hidden.states=TRUE)
#exclude dual transitions between both the observed trait and the hidden trait (e.g., q0A<->q1B)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
trans.rates.hisse.test <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.hisse.test <- ParDrop(trans.rates.hisse.test, c(3,5,8,9,10,12))
trans.rates.hisse.red <- trans.rates.hisse
trans.rates.hisse.red.test <- trans.rates.hisse.test
trans.rates.hisse.red[!is.na(trans.rates.hisse.red) & !trans.rates.hisse.red == 0] = 1
trans.rates.hisse.red.test[!is.na(trans.rates.hisse.red.test) & !trans.rates.hisse.red.test == 0] = 1
##bisse model without hidden states
trans.rates.bisse <- TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse.red <- trans.rates.bisse
trans.rates.bisse.red[!is.na(trans.rates.bisse.red)] = 1
###hisse model with irreversible states
trans.rates.hisse.irrev <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.hisse.irrev <- ParDrop(trans.rates.hisse.irrev, c(1,3,5,8,9,10))
###hisse model with three trans rates, one for 0 -> 1, one for 1 -> 0, and one for transitions among hidden states####
### the following is straight from the HiSSE Vignette:
trans.rates.nodual.threerates <- trans.rates.hisse.red
# Set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving 
# you an estimate of the rate by which shifts in diversification occur):
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3


##Fit the 24 models described in Beaulieu and O'Meara 2016, plus 6 additional models used by Harrington & Reeder (2016)
##The example we use here is topology A and trait=Sting; for running the remaining analyses, change the phy, data and f arguments

hisse.fit1 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse, sann=FALSE)	
hisse.fit2 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse, sann=FALSE)
hisse.fit3 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse.red, sann=FALSE)	
hisse.fit4 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse.red, sann=FALSE)
hisse.fit5 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit6 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit7 <- hisse.null4.old(TopA, Sting, f=sampfrac_s, sann=FALSE)
hisse.fit8 <- hisse.null4.old(TopA, Sting, f=sampfrac_s, eps.anc=rep(1,8), sann=FALSE)
hisse.fit9 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit10 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit11 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit12 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit13 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit14 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit15 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit16 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)	
hisse.fit17 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit18 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)	
hisse.fit19 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit20 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)	
hisse.fit21 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit22 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)	
hisse.fit23 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit24 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)	
hisse_full <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse, sann=FALSE)
hisse_full_irrev <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.irrev, sann=FALSE)
CID2_3rate <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual.threerates, sann=FALSE)
CID2_8rate <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse, sann=FALSE)
CID4_3rate <- hisse.null4.old(TopA, Sting, f=sampfrac_s, trans.type="three.rate", sann=FALSE)
CID4_9rate <- hisse.null4.mod.9.rates(TopA, Sting, f=sampfrac_s, trans.type="All.no.dual", sann=FALSE)


##show the results
hisse.fit1
hisse.fit2
hisse.fit3
hisse.fit4
hisse.fit5
hisse.fit6
hisse.fit7
hisse.fit8
hisse.fit9
hisse.fit10
hisse.fit11
hisse.fit12
hisse.fit13
hisse.fit14
hisse.fit15
hisse.fit16
hisse.fit17
hisse.fit18
hisse.fit19
hisse.fit20
hisse.fit21
hisse.fit22
hisse.fit23
hisse.fit24
hisse_full
hisse_full_irrev
CID2_3rate
CID2_8rate
CID4_3rate
CID4_9rate

###Reconstructing ancestral states for best model ###
CID4_9rate.recon <- MarginRecon.old(TopA, Sting, f = sampfrac_s, four.state.null=T, pars=CID4_9rate$solution, hidden.states=CID4_9rate$hidden.states, AIC=CID4_9rate$AIC, n.cores=2)

###PLOT ancestral states###
recon_plot <- plot.hisse.states(CID4_9rate.recon, rate.param = "net.div", type = "fan", show.tip.label = FALSE, label.offset = 2, fsize = 0.2, legend = "tips", legend.cex = 0.6, edge.width = 2, width.factor=0.4, rate.colors = c("darkcyan", "goldenrod"), state.colors = c("grey", "black"))
