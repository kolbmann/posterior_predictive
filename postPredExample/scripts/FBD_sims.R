library("FossilSim", lib.loc="~/R/win-library/3.4")
library("TreeSim", lib.loc="~/R/win-library/3.4")
data <- read.delim("../TH_rates.txt")
sample <- data[sample(1:nrow(data), 1, replace=FALSE),]

age = sample$origin_time
numbsim = 2
lambda = sample$speciation_rate.1.014369232
mu = sample$extinction_rate.0.862983386
psi = sample$psi
n = 756

# wrong plot
age_tree <- sim.bd.age(age, numbsim, lambda, mu)
plot(age_tree[[1]])

# THIS ONE IS THE PHYLOGENY
taxa_tree <- sim.bd.taxa(n, numbsim, lambda, mu, complete = TRUE, stochsampling = TRUE)
plot(taxa_tree[[1]])

# doesn't finish running
aged_taxa_tree <- sim.bd.taxa.age(n, numbsim, lambda, mu,age = age)
plot(aged_taxa_tree[[1]])

# In if ((lambda/(lambda + mu)) > specevent) { :
# the condition has length > 1 and only the first element will be used
age_FBD_tree <- sim.fbd.age(age, numbsim, lambda, mu, psi, mrca = FALSE, K = 0)
plot(age_FBD_tree[[1]])

# Error in if (any(rate < 0)) stop("Rates must be positive numbers") : 
# missing value where TRUE/FALSE needed
taxa_FBD_tree <- sim.fbd.taxa(n, numbsim, lambda, mu, psi)
plot(taxa_FBD_tree[[1]])
