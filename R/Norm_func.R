# Benchmark Normalization Function Most used: global mean, global median, quantile, lowest CV;

## loess and vsn can be added easily The function also outputs the name of the genes selected by
## lowest CV with $endo

##mean and median uses all data as is.
##quantile uses imputed values NA to largest.  (no Ct or missing criterias used)
##lowest cv and miranorm both use a screened dataset based on Ct and missing criterias. (no imputation)


norm.benchmark = function(data, ct = 25, missing = 0, gene.clust) {

    s.dat <- data[, c("Sample", "Gene", "Ct")]
    colnames(s.dat) <- c("Sample", "Gene", "value")
    dat.cast <- reshape2::acast(s.dat, Gene ~ Sample)

    # Global Median
    normal.med <- normal.gmedian(dat.cast)
    sample.str <- strsplit(as.character(normal.med$Sample), "[.]")
    normal.med$Trt <- unlist(lapply(sample.str, function(x) {x[1]}))
    normal.med$Sample <- unlist(lapply(sample.str, function(x) {x[2]}))

    # Global Mean
    normal.mean <- normal.gmean(dat.cast)
    sample.str <- strsplit(as.character(normal.mean$Sample), "[.]")
    normal.mean$Trt <- unlist(lapply(sample.str, function(x) {x[1]}))
    normal.mean$Sample <- unlist(lapply(sample.str, function(x) {x[2]}))


    # Quantile Normalization
    imp.dat <- impmat(dat.cast)
    normal.qnt <- normal.quant(imp.dat)
    sample.str <- strsplit(as.character(normal.qnt$Sample), "[.]")
    normal.qnt$Trt <- unlist(lapply(sample.str, function(x) {x[1]}))
    normal.qnt$Sample <- unlist(lapply(sample.str, function(x) {x[2]}))


    # Lowest CV normalization Coefficient of variation
    dat <- scr(dat.cast, threshold = ct, missing.percent.thresh = missing)
    dat.sd <- apply(dat, 1, sd, na.rm = TRUE)
    dat.mean <- apply(dat, 1, mean, na.rm = TRUE)
    cv <- dat.sd/dat.mean
    cv.rank <- rank(cv)
    lcv <- names(sort(cv.rank)[1:3])
    normal.hkg <- normal.selected(dat.cast, lcv)
    sample.str <- strsplit(as.character(normal.hkg$Sample), "[.]")
    normal.hkg$Trt <- unlist(lapply(sample.str, function(x) {x[1]}))
    normal.hkg$Sample <- unlist(lapply(sample.str, function(x) {x[2]}))



    # Supplementary function for norm.cluster using gene cluster from miranorm
    dat <- scr(dat.cast, threshold = ct, missing.percent.thresh = missing)
    u.gene <- rownames(dat)
    normal.clust <- normal.selected(dat.cast, gene.clust)
    sample.str <- strsplit(as.character(normal.clust$Sample), "[.]")
    normal.clust$Trt <- unlist(lapply(sample.str, function(x) {x[1]}))
    normal.clust$Sample <- unlist(lapply(sample.str, function(x) {x[2]}))

    # Loess Normalization imp.normal.data = normalize.loess(imp.liver.dat) normal.loess =
    # integ(liver.dat, imp.normal.data)


    ## Variation Stablization fit = vsn2(imp.liver.dat) imp.normal.data0 = predict(fit, newdata =
    ## imp.liver.dat) imp.normal.data = imp.normal.data0 + median(liver.dat$Ct) -
    ## median(imp.normal.data0) normal.vsn = integ(liver.dat, imp.normal.data)

    return(list(normal.mean = normal.mean, normal.median = normal.med, normal.qnt = normal.qnt, normal.hkg = normal.hkg,
        lcv = lcv, normal.clust = normal.clust, u.gene = u.gene))  #, normal.loe = normal.loess, normal.vsn = normal.vsn))
}

normal.gmedian <- function(dat.cast) {
    normal.factor <- apply(dat.cast, 2, median, na.rm = TRUE)
    normal.dat <- sweep(dat.cast, 2, normal.factor, "-")
    dat <- reshape2::melt(normal.dat)
    colnames(dat) <- c("Gene", "Sample", "Ct")
    return(dat)
}


normal.gmean <- function(dat.cast) {
    normal.factor <- apply(dat.cast, 2, mean, na.rm = TRUE)
    normal.dat <- sweep(dat.cast, 2, normal.factor, "-")
    dat <- reshape2::melt(normal.dat)
    colnames(dat) <- c("Gene", "Sample", "Ct")
    return(dat)
}


normal.quant = function(imp.dat) {
    r.names = rownames(imp.dat)
    dmat = order = NULL
    for (i in 1:ncol(imp.dat)) {
        dat = imp.dat[, i]
        o = order(dat)
        odat = dat[o]
        dmat = cbind(dmat, odat)
        order = cbind(order, o)
    }
    for (i in 1:nrow(dmat)) {
        row.mean = mean(dmat[i, ])
        dmat[i, ] = rep(row.mean, ncol(dmat))
    }
    for (j in 1:ncol(dmat)) {
        o = order[, j]
        dmat[, j][o] = dmat[, j]
    }
    colnames(dmat) = colnames(imp.dat)
    rownames(dmat) = r.names
    dat <- reshape2::melt(dmat)
    colnames(dat) <- c("Gene", "Sample", "Ct")
    return(dat)
}


normal.selected <- function(dat.cast, endo) {
    dat.selected <- dat.cast[rownames(dat.cast) %in% endo, ]
    normal.factor <- apply(dat.selected, 2, mean, na.rm = TRUE)
    normal.dat <- sweep(dat.cast, 2, normal.factor, "-")
    dat <- reshape2::melt(normal.dat)
    colnames(dat) <- c("Gene", "Sample", "Ct")
    return(dat)
}

