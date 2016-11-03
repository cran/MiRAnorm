#' Replacing values above ct to NA
#'
#'\code{scr} screens the dataset in wide format according to Ct threshold and missing.percent.thresh (for each gene, the percentage of samples is missing).
#' sets all values above Ct threshold to NA and then removes all rows of genes with number of missing above missing.percent.thresh.
#'
#' @param dat the dataset with rows = Genes and columns = Samples.
#' @param threshold Ct value at or above which all values are set to NA.
#' @param missing.percent.thresh a number between 0 to 100 denoting the percentage of allowable missing.  If the percentage of missing is greater than this value, the row (gene) is not retained.
#'
#' @export


## Screening the dataset in wide format according to the threshold (Ct threshold) and missing.percent.thresh (for each gene, the percentage of samples is missing)


scr = function(dat, threshold, missing.percent.thresh){

    dat[dat>= threshold] <- NA
    rmv.ind = apply(dat, 1, function(x) sum(is.na(x))/length(x) * 100<=missing.percent.thresh)
    dat.new = dat[rmv.ind, ]
    return(data = dat.new)

}

