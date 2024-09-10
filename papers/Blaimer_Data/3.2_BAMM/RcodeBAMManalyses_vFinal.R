This code is for setting up and running diversification analysis in BAMMtools and BAMM, as published in Blaimer et al. XXXX

##Note: this code uses Topology C-1, but input and BAMM results files for Topology A-0 are provided.
HymTopC <- read.tree(file="HYMtopC_bamm_hisse.phy")
setBAMMpriors (HymTopC, total.taxa = 152691, outfile = "HymCpriors.txt", Nmax=1000)

###modify .config with the above returned priors and run BAMM analyses###

###Analyze BAMM output

########Evaluate ESS values/convergence#########
################################################

mcmcout <- read.csv("HYMtopC_bamm_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)




###EVALUATING NUMBER of SHIFTS and best model
##############################################

##Identify number of shifts
edata <- getEventData(HymTopC, eventdata = "HYMtopC_bamm_event_data.txt", burnin=0.1)
shift_probs <- summary(edata)
computeBayesFactors(mcmcout, expectedNumberOfShifts=1, burnin=0.1)
##show results:
shift_probs

######Plotting results#######
###########################

##plot a mean phylorate plot
q <- plot.bammdata(edata, legend=TRUE, lwd=2, method="polar", pal="BrBG", breaksmethod='jenks')
##plot a mean phylorate plot depicting net diversification rates
q2 <- plot.bammdata(edata, legend=TRUE, lwd=2, method="polar", pal="BrBG", breaksmethod='jenks', spex = "netdiv")

##plot histogram of rates
ratesHistogram(q, plotBrks = TRUE, xlab = 'speciation rates')
ratesHistogram(q2, plotBrks = TRUE, xlab = 'net diversification rates')

##plot the overall best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(best, lwd=2, method='polar', legend=TRUE, pal="BrBG", breaksmethod='jenks')
addBAMMshifts(best, method='polar', col='white', bg='black', cex=2)
##plot the overall best shift configuration depicting net diversification rates
plot.bammdata(best, lwd=2, method='polar', legend=TRUE, pal="BrBG", breaksmethod='jenks', spex = "netdiv")
addBAMMshifts(best, method='polar', col='white', bg='black', cex=2)

##summarize and plot credible set of rate shifts
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
summary(css)
plot.credibleshiftset(css, lwd=1.5, pal="BrBG")##compute the marginal shift probabilities

####################################################################################################
#####Other ways to analyze and display rate variation and summarize probabilities of rate shifts####
####################################################################################################

##Plot marginal probability tree, where the branch lengths are scaled by the probability that they contain a shift event
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex=0.5)

#Compute and plot a mean branch tree, where branch lenghts are converted to the mean of the marginal distribution of evolutionary rates on each branch
meanbranchtree_sp <- getMeanBranchLengthTree(edata)
meanbranchtree_ndr <- getMeanBranchLengthTree(edata, rate = "ndr") #to depict net diversification rate instead of speciation rate
plot(meanbranchtree_sp$phy, lwd=2, cex=0.2)
write.tree(meanbranchtree_sp$phy, file="meanbranchtree_TopA.tre")

##Macroevolutionary cohort analysis displays the pairwise probability that any two species share a common macroevolutionary rate dynamic
cmat <- getCohortMatrix(edata)
cohorts(cmat, edata, lwd=1, use.plot.bammdata=TRUE)

##Plot the cumulative shift probability tree, showing the cumulative probability, on each branch, that a shift occurred somewhere between the focal branch and the root of the tree
cst <- cumulativeShiftProbsTree(edata)
plot.phylo(cst, cex=0.2)

cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('gray48', length(HymTopAu$edge.length))
is_highprobshift <- cst$edge.length >= 0.97
is_highestprobshift <- cst$edge.length >= 1.00
edgecols[ is_highprobshift ] <- "darkcyan"
edgecols[ is_highestprobshift ] <- "goldenrod3"
plot.phylo(HymTopAu, edge.color = edgecols, type="fan",show.tip.label = F)


##################################################
###compute mean speciation and extinction rates for focal clades with indicated rate shifts### 

#Hymenoptera#
allrates <- getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))
mean(allrates$mu)
quantile(allrates$mu, c(0.05, 0.95))

#Vespina#
vespinarates <- getCladeRates(edata, node=771)
mean(vespinarates$lambda)
quantile(vespinarates$lambda, c(0.05, 0.95))
mean(vespinarates$mu)
quantile(vespinarates$mu, c(0.05, 0.95))
nonvespinarate <- getCladeRates(edata, node = 771, nodetype = "exclude")
mean(nonvespinarate$lambda)
quantile(nonvespinarate$lambda, c(0.05, 0.95))
mean(nonvespinarate$mu)
quantile(nonvespinarate$mu, c(0.05, 0.95))


#Xiphydriidae + Siricidae#
woodwasprates <- getCladeRates(edata, node=1509)
mean(woodwasprates$lambda)
quantile(woodwasprates$lambda, c(0.05, 0.95))
mean(woodwasprates$mu)
quantile(woodwasprates$mu, c(0.05, 0.95))
nonwoodwasprate <- getCladeRates(edata, node = 1509, nodetype = "exclude")
mean(nonwoodwasprate$lambda)
quantile(nonwoodwasprate$lambda, c(0.05, 0.95))
mean(nonwoodwasprate$mu)
quantile(nonwoodwasprate$mu, c(0.05, 0.95))

#Pamphiliidae + Xyelidae#
PampXyerates <- getCladeRates(edata, node=1511)
mean(PampXyerates$lambda)
quantile(PampXyerates$lambda, c(0.05, 0.95))
mean(PampXyerates$mu)
quantile(PampXyerates$mu, c(0.05, 0.95))
nonPampXyerate <- getCladeRates(edata, node = 1509, nodetype = "exclude")
mean(nonPampXyerate$lambda)
quantile(nonPampXyerate$lambda, c(0.05, 0.95))
mean(nonPampXyerate$mu)
quantile(nonPampXyerate$mu, c(0.05, 0.95))

#Aculeata excl. Chrysidoidea#
aculeaterates <- getCladeRates(edata, node=1060)
mean(aculeaterates$lambda)
quantile(aculeaterates$lambda, c(0.05, 0.95))
mean(aculeaterates$mu)
quantile(aculeaterates$mu, c(0.05, 0.95))
nonaculeaterate <- getCladeRates(edata, node = 1060, nodetype = "exclude")
mean(nonaculeaterate$lambda)
quantile(nonaculeaterate$lambda, c(0.05, 0.95))
mean(nonaculeaterate$mu)
quantile(nonaculeaterate$mu, c(0.05, 0.95))

#Bees excl. Melittidae#
beesrates <- getCladeRates(edata, node=1069)
mean(beesrates$lambda)
quantile(beesrates$lambda, c(0.05, 0.95))
mean(beesrates$mu)
quantile(beesrates$mu, c(0.05, 0.95))
nonbeerate <- getCladeRates(edata, node = 1069, nodetype = "exclude")
mean(nonbeerate$lambda)
quantile(nonbeerate$lambda, c(0.05, 0.95))
mean(nonbeerate$mu)
quantile(nonbeerate$mu, c(0.05, 0.95))

#Cynipidae s.s. + Figitidae s.l.#
cynifigirates <- getCladeRates(edata, node=983)
mean(cynifigirates$lambda)
quantile(cynifigirates$lambda, c(0.05, 0.95))
mean(cynifigirates$mu)
quantile(cynifigirates$mu, c(0.05, 0.95))
noncynifigirate <- getCladeRates(edata, node = 983, nodetype = "exclude")
mean(noncynifigirate$lambda)
quantile(noncynifigirate$lambda, c(0.05, 0.95))
mean(noncynifigirate$mu)
quantile(noncynifigirate$mu, c(0.05, 0.95))

#Cynipoidea#
cynipoidrates <- getCladeRates(edata, node=980)
mean(cynipoidrates$lambda)
quantile(cynipoidrates$lambda, c(0.05, 0.95))
mean(cynipoidrates$mu)
quantile(cynipoidrates$mu, c(0.05, 0.95))
noncynipoidrate <- getCladeRates(edata, node = 980, nodetype = "exclude")
mean(noncynipoidrate$lambda)
quantile(noncynipoidrate$lambda, c(0.05, 0.95))
mean(noncynipoidrate$mu)
quantile(noncynipoidrate$mu, c(0.05, 0.95))

#Ichneumonidae, internal#
Ichrates <- getCladeRates(edata, node=1281)
mean(Ichrates$lambda)
quantile(Ichrates$lambda, c(0.05, 0.95))
mean(Ichrates$mu)
quantile(Ichrates$mu, c(0.05, 0.95))
nonIchrate <- getCladeRates(edata, node = 1281, nodetype = "exclude")
mean(nonIchrate$lambda)
quantile(nonIchrate$lambda, c(0.05, 0.95))
mean(nonIchrate$mu)
quantile(nonIchrate$mu, c(0.05, 0.95))

#Ichneumonoidea#
Ichnoidrates <- getCladeRates(edata, node=1212)
mean(Ichnoidrates$lambda)
quantile(Ichnoidrates$lambda, c(0.05, 0.95))
mean(Ichnoidrates$mu)
quantile(Ichnoidrates$mu, c(0.05, 0.95))
nonIchnoidrate <- getCladeRates(edata, node = 1212, nodetype = "exclude")
mean(nonIchnoidrate$lambda)
quantile(nonIchnoidrate$lambda, c(0.05, 0.95))
mean(nonIchnoidrate$mu)
quantile(nonIchnoidrate$mu, c(0.05, 0.95))

#Tenthredinoidea#
Tentrates <- getCladeRates(edata, node=1514)
mean(Tentrates$lambda)
quantile(Tentrates$lambda, c(0.05, 0.95))
mean(Tentrates$mu)
quantile(Tentrates$mu, c(0.05, 0.95))
nonTentrate <- getCladeRates(edata, node = 1514, nodetype = "exclude")
mean(nonTentrate$lambda)
quantile(nonTentrate$lambda, c(0.05, 0.95))
mean(nonTentrate$mu)
quantile(nonTentrate$mu, c(0.05, 0.95))

#Apocrita, excl. Ichneumonoidea#
Apoexrates <- getCladeRates(edata, node=773)
mean(Apoexrates$lambda)
quantile(Apoexrates$lambda, c(0.05, 0.95))
mean(Apoexrates$mu)
quantile(Apoexrates$mu, c(0.05, 0.95))
nonApoexrate <- getCladeRates(edata, node = 773, nodetype = "exclude")
mean(nonApoexrate$lambda)
quantile(nonApoexrate$lambda, c(0.05, 0.95))
mean(nonApoexrate$mu)
quantile(nonApoexrate$mu, c(0.05, 0.95))

#Apocrita#
Aporates <- getCladeRates(edata, node=772)
mean(Aporates$lambda)
quantile(Aporates$lambda, c(0.05, 0.95))
mean(Aporates$mu)
quantile(Aporates$mu, c(0.05, 0.95))
nonAporate <- getCladeRates(edata, node = 772, nodetype = "exclude")
mean(nonAporate$lambda)
quantile(nonAporate$lambda, c(0.05, 0.95))
mean(nonAporate$mu)
quantile(nonAporate$mu, c(0.05, 0.95))

#Cynipidae ss#
cynirates <- getCladeRates(edata, node=984)
mean(cynirates$lambda)
quantile(cynirates$lambda, c(0.05, 0.95))
mean(cynirates$mu)
quantile(cynirates$mu, c(0.05, 0.95))
noncynirate <- getCladeRates(edata, node = 984, nodetype = "exclude")
mean(noncynirate$lambda)
quantile(noncynirate$lambda, c(0.05, 0.95))
mean(noncynirate$mu)
quantile(noncynirate$mu, c(0.05, 0.95))

#Eurytomidae, internal topC#
euryrates <- getCladeRates(edata, node=824)
mean(euryrates$lambda)
quantile(euryrates$lambda, c(0.05, 0.95))
mean(euryrates$mu)
quantile(euryrates$mu, c(0.05, 0.95))
noneuryrate <- getCladeRates(edata, node = 824, nodetype = "exclude")
mean(noneuryrate$lambda)
quantile(noneuryrate$lambda, c(0.05, 0.95))
mean(noneuryrate$mu)
quantile(noneuryrate$mu, c(0.05, 0.95))

#Eurytomidae, internal topA#
eury2rates <- getCladeRates(edata, node=813)
mean(eury2rates$lambda)
quantile(eury2rates$lambda, c(0.05, 0.95))
mean(eury2rates$mu)
quantile(eury2rates$mu, c(0.05, 0.95))
noneury2rate <- getCladeRates(edata, node = 813, nodetype = "exclude")
mean(noneury2rate$lambda)
quantile(noneury2rate$lambda, c(0.05, 0.95))
mean(noneury2rate$mu)
quantile(noneury2rate$mu, c(0.05, 0.95))

#Formicidae#
formirates <- getCladeRates(edata, node=1128)
mean(formirates$lambda)
quantile(formirates$lambda, c(0.05, 0.95))
mean(formirates$mu)
quantile(formirates$mu, c(0.05, 0.95))
nonformirate <- getCladeRates(edata, node = 1128, nodetype = "exclude")
mean(nonformirate$lambda)
quantile(nonformirate$lambda, c(0.05, 0.95))
mean(nonformirate$mu)
quantile(nonformirate$mu, c(0.05, 0.95))

