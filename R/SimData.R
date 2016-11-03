#' Simulating a RT-PCR miRNA dataset.
#'
#' \code{simData} simulates a RT-PCR miRNA dataset with user defined levels of variability and treatment effect size.
#'
#' @param n.trt Number of simulated treatment samples
#' @param n.ctrl Number of simulated control samples
#' @param n.gene Number of simulated genes in the panel
#' @param n.big.effect Number of genes with large treatment effect
#' @param n.small.effect Number of genes with small treatment effect
#' @param mean.big.effect Average effect size for a "large" treatment effect
#' @param mean.small.effect Average effect size for a "small" treatment effect
#' @param sigma.gene sd of gene to gene effect sizes for large and small treatment effects.
#' @param sigma.error a vector of length 2 for 2 different measurement error sizes (sd).
#' @param n.err number of genes with sigma.error[2].  The rest (n.gene - n.err) have sigma.error[1] measurement error.
#' @param mean.sample the unadjusted mean of the samples.  Can generally be left as default
#' @param sigma.sample the unadjusted sample to sample sd.  Can generallyb e left as default.
#'
#' @export


###Simulation####

###some input criterion####


simData<- function(n.trt=50, n.ctrl=50, n.gene=96, sigma.error, n.err = 10, mean.sample=0.6, sigma.sample = 0.6, sigma.gene = 0,n.big.effect = 5, n.small.effect = 15, mean.big.effect = 5,mean.small.effect = 2)
{
  n.no.effect <- n.gene - n.big.effect - n.small.effect; ###remainder genes with no effect###
  ###size of effects in trt grp###
  size.big.effect <- mean.big.effect + rnorm(n.big.effect, sd=sigma.gene);
  size.small.effect <- mean.small.effect + rnorm(n.small.effect, sd=sigma.gene);
  #size.big.effect <- sample(c(1,-1), size=n.big.effect, replace=TRUE, prob = c(0.5, 0.5)) * (mean.big.effect + rnorm(n.big.effect, sd=sigma.gene));
  #size.small.effect <- sample(c(1,-1), size=n.small.effect, replace=TRUE, prob = c(0.5, 0.5)) * (mean.small.effect + rnorm(n.small.effect, sd=sigma.gene));
  ###true underlying latent values###
  true.gene.ctrl <- runif(n.gene, min=10, max=20);  ##mean of each gene in the control group###
  true.gene.trt <-  true.gene.ctrl + c(size.big.effect, size.small.effect, rep(0, n.no.effect))

  ##sample.level differences

  sample.eff.trt <- rnorm(n = n.trt, sd = sigma.sample) + sample(c(0, 1), n.trt, replace = TRUE, prob = c(0.7, 0.3))* mean.sample;

  sample.eff.ctrl <- rnorm(n = n.ctrl, sd = sigma.sample) + sample(c(0, 1), n.ctrl, replace = TRUE, prob = c(0.7, 0.3))* mean.sample;

  ##different level of measurement error
  ## measurement error is considered in a vector format, deciding how much of data is with the specified measurement error
  measure.err.t = measure.err.c = NULL
  #n.sigma =  n.gene %/% length(sigma.error)
  for (i in 1:n.trt){
      measure.err.t <- c(measure.err.t, rnorm(n = n.gene - n.err, sd = sigma.error[1]), rnorm(n = n.err, sd = sigma.error[2]))
  }
  for (i in 1:n.ctrl){
      measure.err.c <- c(measure.err.c, rnorm(n = n.gene - n.err, sd = sigma.error[1]), rnorm(n = n.err, sd = sigma.error[2]))
  }

  ###simulate ctrl grp data###

  dat.ctrl <- data.frame(Ct = rep(true.gene.ctrl, times=n.ctrl) + rep(sample.eff.ctrl, each=n.gene) + measure.err.c,
                         Gene = rep(paste("Gene", 1:n.gene, sep="."), times=n.ctrl),
                         Sample = rep(1:n.ctrl, each=n.gene),
                         Group = rep("Ctrl", n.ctrl*n.gene))

  dat.trt <- data.frame(Ct = rep(true.gene.trt, times=n.trt) + rep(sample.eff.trt, each=n.gene) + measure.err.t,
                        Gene = rep(paste("Gene", 1:n.gene, sep="."), times=n.trt),
                        Sample = rep((n.ctrl+1):(n.ctrl+n.trt), each=n.gene),
                        Group = rep("Trt", n.trt*n.gene))

  sim.dat <- rbind(dat.ctrl, dat.trt)
  return(list(sim = sim.dat, trt.eff= true.gene.trt[1: (n.small.effect+n.big.effect+n.no.effect)]-true.gene.ctrl[1: (n.small.effect+n.big.effect+n.no.effect)]))
}

