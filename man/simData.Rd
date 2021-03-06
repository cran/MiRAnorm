% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimData.R
\name{simData}
\alias{simData}
\title{Simulating a RT-PCR miRNA dataset.}
\usage{
simData(n.trt = 50, n.ctrl = 50, n.gene = 96, sigma.error, n.err = 10,
  mean.sample = 0.6, sigma.sample = 0.6, sigma.gene = 0,
  n.big.effect = 5, n.small.effect = 15, mean.big.effect = 5,
  mean.small.effect = 2)
}
\arguments{
\item{n.trt}{Number of simulated treatment samples}

\item{n.ctrl}{Number of simulated control samples}

\item{n.gene}{Number of simulated genes in the panel}

\item{sigma.error}{a vector of length 2 for 2 different measurement error sizes (sd).}

\item{n.err}{number of genes with sigma.error[2].  The rest (n.gene - n.err) have sigma.error[1] measurement error.}

\item{mean.sample}{the unadjusted mean of the samples.  Can generally be left as default}

\item{sigma.sample}{the unadjusted sample to sample sd.  Can generallyb e left as default.}

\item{sigma.gene}{sd of gene to gene effect sizes for large and small treatment effects.}

\item{n.big.effect}{Number of genes with large treatment effect}

\item{n.small.effect}{Number of genes with small treatment effect}

\item{mean.big.effect}{Average effect size for a "large" treatment effect}

\item{mean.small.effect}{Average effect size for a "small" treatment effect}
}
\description{
\code{simData} simulates a RT-PCR miRNA dataset with user defined levels of variability and treatment effect size.
}

