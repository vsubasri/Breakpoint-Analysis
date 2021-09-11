library(doParallel)
library(plyr)
library(dotCall64)

registerDoParallel(cores = 3)
genome = 31371612 #3137161264
N = 10

setwd("~/Medical Biophysics 2017:2018/Rotation Project 3")
mutation_rates <- read.delim("mutation_rates.tsv")

run_simulations <- function(N, rate, sample, pre) {
   mylist <- list()
   locations <- list()
   post <- pre
   for(i in 1:N) {
     message(paste("Simulation: ", i))
     sim <- rbinom(n = genome, size = 1, prob = rate)
     indices <- which(sim %in% c(1))
     for(ind in indices) {
        post[ind] = sample(c('A','T','G','C'), size=1, replace=TRUE)
     }
     combined <- paste0(pre, post)
     mylist[[i]] <- as.data.frame.list(table(combined))
     locations[[i]] <- as.data.frame(t(indices))
   }
   df <- do.call("rbind.fill",mylist)
   loc <- do.call("rbind.fill",locations)
   write.table(df, file = paste0("simulations_",sample,".tsv"), sep='\t', quote=FALSE)
   write.table(loc, file = paste0("simlocations_",sample,".tsv"), sep='\t', quote=FALSE)
   #save(df, file = paste0("simulations_",sample,".rda"), compress = "xz")
}

run_all <- function(N) {
  for(i in 1:dim(mutation_rates)[1]) {
    rate = mutation_rates$mutation.rate[i]
    message(mutation_rates$sample[i])
    sample <- mutation_rates$sample[i]
    pre <- sample(c('A','T','G','C'), size=genome, replace=TRUE)
    run_simulations(N, rate, sample, pre)
    break
  }
}

run_all(N)


