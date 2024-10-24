extract_last_number <- function(state) {
  if (grepl("\\(.*\\)", state)) {  # Check if parentheses are present
    # Extract the numbers and return the last one
    numbers <- unlist(regmatches(state, gregexpr("[0-9]+", state)))
    return(as.numeric(numbers[length(numbers)]))  # Return the last number
  }
  return(as.numeric(state))  # If no parentheses, return the state as numeric
}


library(diversitree)
library(phytools)

setwd("..")


data =  read.csv("data/final_MuSSE_chars_4_5.csv")

tree = read.tree("data/pruned_dasilva_tree.tre")

tree = force.ultrametric(tree)

# we have 4 caracter states, what does the MuSSE arguement list look like

diversitree:::default.argnames.musse(4)


data[,2] = gsub("\n    ", "", data[,2])


data[,2] = unlist(lapply( 1:length(data$State), function(i) extract_last_number( data$State[i])))

tip.states = data[,2]

names(tip.states) = data$Species

pars <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
set.seed(2)
phy <- tree.musse(pars, 20, x0 = 1)


#lik.b <- make.bisse(phy, phy$tip.state - 1)
lik.m <- make.musse(tree, tip.states, sampling.f=c(.01, .01, .01, .01),  4)


fit <- find.mle(lik.m, pars, method="subplex")

fit$convergence

named_vector = fit$par

# Convert the named vector into a data frame for boxplotting
df <- data.frame(
  names = names(named_vector),
  values = as.numeric(named_vector)
)

# Create a boxplot
boxplot(values ~ names, data = df, 
        main = "Boxplot of Named Vector Elements", 
        xlab = "Names", 
        ylab = "Values",
        col = "lightblue", 
        border = "black")



all.equal(lik.m(pars), lik.m(pars), tolerance = 1e-07)





diversitree:::default.argnames.musse(4)


pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
