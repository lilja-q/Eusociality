# This script fits the MiSSE model, a completely trait-free version of a HiSSE model
# This code is taken from the MiSSE vignette by Jeremy Beaulieu within the HiSSE R package.

#Load relevant packages
library(hisse)

#load tree file and data 

TopA <- read.tree(file="HYMtopA_bamm_hisse.phy")
TopC <- read.tree(file="HYMtopC_bamm_hisse.phy")
f=0.005 #sampling frequency

##fit the models
TopAmisse <- MiSSEGreedy(TopA, f=0.005, possible.combos = generateMiSSEGreedyCombinations(), stop.deltaAICc=10, save.file=NULL, n.cores=6, chunk.size=10, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=TRUE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0)
TopCmisse <- MiSSEGreedy(TopC, f=0.005, possible.combos = generateMiSSEGreedyCombinations(), stop.deltaAICc=10, save.file=NULL, n.cores=6, chunk.size=10, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=TRUE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0)

##Summarize results for all models## Note: This computationally expensive and not absolutely necessary to move on; one can simply extract and compare model scores from the MiSSEGreedy objects.
TopAmisse_states <- SummarizeMiSSEGreedy(TopAmisse, min.weight=0.01, n.cores=6, recon=TRUE)
TopCmisse_states <- SummarizeMiSSEGreedy(TopCmisse, min.weight=0.01, n.cores=6, recon=TRUE)

##Reconstruct the states using best model (replace "#best_model" with number of highest-scoring model)
TopAmisse_states_best <- MarginReconMiSSE(TopA, f=0.005, TopAmisse[[#best_model]]$solution, TopAmisse[[#best_model]]$hidden.states, fixed.eps=TopAmisse[[#best_model]]$fixed.eps, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, AIC=TopAmisse[[#best_model]]$AIC, get.tips.only=FALSE, verbose=TRUE, n.cores=6, dt.threads=1)
TopCmisse_states_best <- MarginReconMiSSE(TopC, f=0.005, TopCmisse[[#best_model]]$solution, TopCmisse[[#best_model]]$hidden.states, fixed.eps=TopCmisse[[#best_model]]$fixed.eps, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, AIC=TopCmisse[[#best_model]]$AIC, get.tips.only=FALSE, verbose=TRUE, n.cores=6, dt.threads=1)

##PLOT ancestral states
recon_plot_TopAmisse_best <- plot.misse.states(TopAmisse_states_best, rate.param = "net.div", type = "fan", show.tip.label = TRUE, label.offset = 1, fsize = 0.2, legend = "tips", legend.cex = 0.6, edge.width = 2, width.factor=0.4, rate.colors = c("darkcyan", "goldenrod"), state.colors = c("grey", "black"))
recon_plot_TopCmisse_best <- plot.misse.states(TopCmisse_states_best, rate.param = "net.div", type = "fan", show.tip.label = TRUE, label.offset = 1, fsize = 0.2, legend = "tips", legend.cex = 0.6, edge.width = 2, width.factor=0.4, rate.colors = c("darkcyan", "goldenrod"), state.colors = c("grey", "black"))

##Extract tip rates from best model
tip.rates_TopA_best <- GetModelAveRates(TopAmisse_states_best, type = c("tips"))
tip.rates_TopC_best <- GetModelAveRates(TopCmisse_states_best, type = c("tips"))