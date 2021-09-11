
path<- "~/research/breakpoint_analysis/monte_carlo/whole_genome_simulations"
setwd(path)

N = 1000

mutations <- function(file) {
  message("Computing distances for ", file)
  locs <- read.delim(file)
  sorted <- t(apply(locs, 1, sort)) 
  counter <- 0
  count <- list()
  for(i in 1:length(sorted)) {
    message("Simulation ", i)
    d <- c(dist(unlist(sorted[i])))
    ###
    pw_d <- dist(unlist(sorted[[i]]))
    pw_d_mat <- as.matrix(pw_d, labels=TRUE)
    colnames(pw_d_mat) <- rownames(pw_d_mat) <- sorted[[i]]
    which(pw_d_mat<1000)
    ###
    tot <- length(d)
    lt <- length(d[d < N])
    p <- as.double(lt)
    #    p <- as.double(lt/tot)
    inc(counter)
    count[[counter]] <- p
  }
  count<- unlist(count)
  if (any(is.nan(count))) {
    count[is.nan.data.frame(count)] <- -1
  }
  return(count)
}
  
plot_general <- function(count, sample) {
  c <- as.data.frame(count)
  message(max(table(count))+10)
  hist <- ggplot(c, aes(x=c$count)) + geom_histogram(bins=length(unique(count)+5)) + 
    labs(x = "Number of SNVs (x) within 1kb", y="f(x) based on simulations") + 
    scale_x_continuous(breaks=seq(0, max(count)+2, 1)) + 
    scale_y_continuous(breaks=seq(0, max(table(count))*2, 10)) +
    theme_bw() +
    labs(title=paste0("Monte Carlo Simulation of SNVs occuring \n within 1kb - ",sample)) +
    theme(plot.title = element_text(hjust = 0.5))
  hist
  ggsave(filename=sample, plot = hist, device="pdf")
}
  
plot_dens <- function(count) {
  plot(density(count), main="Monte Carlo Simulation Mutation Rate")
}

plot_dens_gapped <- function(count) {
  dens <- density(count)
  par(mfrow=2:1, mar=rep(.5, 4), oma=c(3,3,1,1), cex=0.5)
  plot(dens, ylim=c(0.2, 1), axes=FALSE, main="")
  axis(2); box()
  lines(dens)
  plot(dens, ylim=c(0, 0.2), xlim =c(0,30), main="")
  lines(dens)
}
  
is.nan.data.frame <- function(x) {
      do.call(cbind, lapply(x, is.nan))
}
  
inc <- function(x)
{
  eval.parent(substitute(x <- x + 1))
  
}

compute_all <- function() {
  all <- list()
  files <- list.files(path, pattern=".tsv", all.files=FALSE, full.names=FALSE) 
  for(i in 1:length(files)) {
    all[[i]] <- mutations(files[i])
  }
  df <- as.data.frame(all, col.names=files)
  return(df)
}

all <- compute_all()

