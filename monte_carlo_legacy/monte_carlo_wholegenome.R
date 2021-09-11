library(doParallel)
library(plyr)

registerDoParallel(cores = 3)
genome = 3137161264
N = 1000

setwd("~/research/breakpoint_analysis")
mutation_rates <- read.delim("mutation_rates.tsv")

run_simulations <- function(N, rate, sample) {
  mylist <- list()
  for(i in 1:N) {
    message(paste("Simulation: ", i))
    sim <- rbinom(n = genome, size = 1, prob = rate)
    sim1 <- which(sim[1:1000000000]==1)
    sim2 <- which(sim[1000000000:2000000000]==1)+1000000000
    sim3 <- which(sim[2000000000:3000000000]==1)+2000000000
    sim4 <- which(sim[3000000000:genome]==1)+3000000000
    # for (ind in 1:length(sim)) {
    #   if (sim[ind] == 1) {
    #     indices <- c(indices, ind)
    #   }
    # }
    indices <- c(sim1, sim2, sim3, sim4)
    for(ind in indices) {
      post[ind] = sample(c('A','T','G','C'), size=1, replace=TRUE)
    }
    mylist[[i]] <- as.data.frame(t(indices))
  }
  locs <- do.call("rbind.fill",mylist)
  write.table(locs, file = paste0("simulations_", sample,".tsv"), sep='\t', quote = FALSE)
}

run_all <- function(N) {
  for(i in 1:dim(mutation_rates)[1]) {
    rate = mutation_rates$mutation.rate[i]
    message(mutation_rates$sample[i])
    sample <- mutation_rates$sample[i]
    run_simulations(N, rate, sample)
  }
}

run_all(N)


