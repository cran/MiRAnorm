#Cluster Selection Algorithm
#
#calculate the whole trace of selection probability over different size levels


### example### obj = trace.selection.kappa(dat, id.var='Sample', n.boot=200, max.clust=20,
### min.clust = 4, dis.method = 'Euclidean', hclust.method = 'single', plot = FALSE, selected = 10,
### multicore = 10)

### n.boot is the bootstrap sample min.clust is the smallest cluster size you consider and
### max.clust is the maximum cluster size you consider. They will determine the range of cluster
### size and the range will affect the magnitude of summary statistic but may not affect the order
### of summary statistic across different miRNAs the specified 'selected' here is only relevant to
### the plot whether the curve is colored or grey co.mat keeps the stability statistic for each
### miRNA over the whole cluster size range sentences below can be used to print out the ordered
### miRNAs according to the stability summary statistic co.mat <- obj$co co.mean <- apply(co.mat,
### 1, mean) co.order <- co.mean[order(co.mean, decreasing = TRUE)] kappa.mat keeps the kappa
### coefficient over the whole cluster size range, s.co.mat is the adjusted co.mat by kappa
### coefficient.

trace.selection.kappa <- function(dat, id.var = "Sample", n.boot, ct = 25, missing = 0, dis.method = "Euclidean",
    hclust.method = "single", min.clust = 3, max.clust = 15, plot = TRUE, plot.more = FALSE, selected = 10,
    multicore = 7) {

  ##included for R cmd check of undefined global variables, note: these are actually defined names
  ##since the datasets intended for the package are well defined and structured.
  Col = Gene = Size = Stability = NULL;


    u.gene = unique(dat$Gene)

    co.mat <- matrix(nrow = length(u.gene), ncol = max.clust - min.clust + 1)
    rownames(co.mat) = as.character(u.gene)
    kappa.mat <- matrix(1, ncol = max.clust - min.clust + 1)


    for (k in min.clust:max.clust) {
        obj = cluster.kappa(dat, id.var = id.var, n.boot = n.boot, size = k, ct = ct, missing = missing,
            dis.method = dis.method, hclust.method = hclust.method, multicore = multicore)
        co.mat[, k - min.clust + 1] = obj$co
        kappa.mat[1, k - min.clust + 1] = obj$kappa
    }

    colnames(co.mat) = as.character(seq(min.clust, max.clust, 1))

    s.co.mat = t(co.mat/as.numeric(kappa.mat))

    if (plot == TRUE) {
        pdf("selection.pdf")
        co.mean <- apply(co.mat, 1, mean)
        co.order <- co.mean[order(co.mean, decreasing = TRUE)]
        co.long <- reshape2::melt(co.mat, id.vars = rownames(co.mat), measure.vars = colnames(co.mat), variable.name = "Stability")
        colnames(co.long) = c("Gene", "Size", "Stability")
        co.long$Col = "Non_selected"
        co.long$Col[co.long$Gene %in% names(co.order[1:selected])] = unique(as.character(co.long$Gene[co.long$Gene %in%
            names(co.order[1:selected])]))
        col.vec = c("dark blue", "blue", "dark red", "red", "brown", "purple", "pink", "orange",
            "yellow", "dark green", "green", "gold", "black")
        print(ggplot(co.long, aes(x = Size, y = Stability, group = Gene, col = Col)) + geom_line(alpha = 0.5) +
            labs(title = "Stability") + theme_bw() + scale_colour_manual(values = c(col.vec[1:selected],
            "grey")))
        # t.co.mat = t(co.mat) matplot(seq(min.clust, max.clust, 1), jitter(t.co.mat), xlab =
        # 'cluster.size', ylab = 'stability', type = 'l', lty = 1, col = rainbow(ncol(t.co.mat)))
        # legend('right', colnames(t.co.mat),col=rainbow(ncol(t.co.mat)),fill=rainbow(ncol(t.co.mat)),
        # cex = 0.3, pt.cex = 1)
        dev.off()
    }
    if (plot.more == TRUE) {
        pdf("kappa.pdf")
        plot(seq(min.clust, max.clust, 1), kappa.mat)
        dev.off()

        pdf("kappa_selection.pdf")
        s.co.mat = t(co.mat/as.numeric(kappa.mat))
        matplot(seq(min.clust, max.clust, 1), jitter(s.co.mat), xlab = "cluster.size", ylab = "stability",
            type = "l", lty = 1, col = rainbow(ncol(s.co.mat)))
        legend("right", colnames(s.co.mat), col = rainbow(ncol(s.co.mat)), fill = rainbow(ncol(s.co.mat)),
            cex = 0.3, pt.cex = 1)
        dev.off()
    }


    return(list(co = co.mat, kappa = kappa.mat, s.co = s.co.mat))

}



### function to generate a list of indices for bootstrapped data###
bt.strap <- function(N, x, n.boot) {
    lapply(1:n.boot, FUN = function(interation, x, N) {
        sample(x = x, size = N, replace = TRUE)
    }, N = N, x = x)
}


### function to calculate the norm.cluster for the specific bootstrapped sample###
index.pair <- function(dat, id.var, dis.method, hclust.method, ct, missing, size, bt.ind) {

    mySubsets <- lapply(bt.ind, function(x, currentDF) {
        currentDF[currentDF[, id.var] == x, ]
    }, dat)
    myNames <- lapply(1:length(mySubsets), function(x, currentDat) {
        currentDat[[x]]$unique = x
        return(currentDat[[x]])
    }, mySubsets)

    dat.resample <- do.call("rbind", myNames)
    dat.resample[, id.var] <- paste(dat.resample[, id.var], dat.resample$unique, sep = ".")
    dat.resample <- dat.resample[, !names(dat.resample) %in% "unique"]

    ncls <- norm.cluster(data = dat.resample, dis.method = dis.method, hclust.method = hclust.method,
        ct = ct, missing = missing, size = size)
    ct.mat <- data.frame(names = unique(dat$Gene), ifelse(unique(dat$Gene) %in% ncls$gene.clust,
        1, 0))

    return(list(cut = ncls$cut.group, ct.mat = ct.mat))

}


### function to calculate Kappa coefficient of two sets based on selection performance###
agree.twosets = function(aset1, aset2, p.tot) {
    if (length(aset1) == 0 || length(aset2) == 0 || length(aset1) + length(aset2) == 2 * p.tot)
        ratio = -1 else {
        n11 = length(intersect(aset1, aset2))
        n22 = length(intersect(setdiff(c(1:p.tot), aset1), setdiff(c(1:p.tot), aset2)))
        n12 = length(setdiff(aset1, intersect(aset1, aset2)))
        n21 = length(setdiff(aset2, intersect(aset1, aset2)))
        ratio = ((n11 + n22)/p.tot - ((n11 + n12) * (n11 + n21) + (n12 + n22) * (n21 + n22))/(p.tot *
            p.tot))/(1 - ((n11 + n12) * (n11 + n21) + (n12 + n22) * (n21 + n22))/(p.tot * p.tot))
    }
    return(ratio)
}


### function to calculate how much level of matching between the bootstrap pairs###
cluster.m <- function(m, dat, id.var, size, ct, missing, dis.method, hclust.method, bt.ind1, bt.ind2) {
    co <- matrix(rep(0, length(unique(dat$Gene))), nrow = length(unique(dat$Gene)))
    rownames(co) <- unique(dat$Gene)

    ct.mat1 = index.pair(dat, id.var, dis.method, hclust.method, ct, missing, size, bt.ind1[[m]])$ct.mat
    ct.mat2 = index.pair(dat, id.var, dis.method, hclust.method, ct, missing, size, bt.ind2[[m]])$ct.mat

    subset = intersect(ct.mat1$names, ct.mat2$names)

    aset1 = ct.mat1[ct.mat1$names %in% subset, 2]
    aset2 = ct.mat2[ct.mat2$names %in% subset, 2]

    dummied = (aset1 + aset2)/2

    co[rownames(co) %in% subset] = dummied

    kappa = agree.twosets(which(aset1 != 0), which(aset2 != 0), length(subset))

    return(list(co = co, kappa = kappa))

}


### function to calculate the overall kappa statistic at certain size level (parallel computation
### for paired samples)###
cluster.kappa <- function(dat, id.var = "Sample", n.boot, ct = 25, missing = 0, size, dis.method = "Euclidean",
    hclust.method = "single", multicore = 4) {

    subj.ID <- dat[, id.var]
    bt.ind1 <- bt.strap(N = length(unique(subj.ID)), x = unique(subj.ID), n.boot = n.boot)
    bt.ind2 <- bt.strap(N = length(unique(subj.ID)), x = unique(subj.ID), n.boot = n.boot)

    ##Check if cores can be forked
    if(.Platform$OS.type == "unix") {
      cl <- parallel::makeCluster(multicore, type = "FORK")
    } else {
      cl <- parallel::makeCluster(multicore, type = "PSOCK")
      parallel::clusterExport(cl, varlist = c("index.pair", "norm.cluster", "n.boot", "cluster.m", "dat", "id.var", "size",
                                              "ct", "missing", "dis.method", "hclust.method", "bt.ind1", "bt.ind2"), envir = environment())
    }

    all.boot <- parallel::parLapply(cl, 1:n.boot, cluster.m, dat, id.var, size, ct, missing, dis.method, hclust.method,
        bt.ind1, bt.ind2)

    parallel::stopCluster(cl)

    # for troubleshooting only
    #all.boot <- lapply(1:n.boot, cluster.m, dat, id.var, size, ct , missing, dis.method, hclust.method, bt.ind1, bt.ind2)
    kappa = 0

    co <- matrix(rep(0, length(unique(dat$Gene))), nrow = length(unique(dat$Gene)))
    rownames(co) <- unique(dat$Gene)

    for (m in 1:n.boot) {

        co = co + all.boot[[m]]$co
        kappa = kappa + all.boot[[m]]$kappa

    }

    kappa.prob = as.numeric(kappa/n.boot)
    # print(kappa.prob)
    co.prob = as.numeric(co/n.boot)
    # print(co.prob)
    return(list(co = co.prob, kappa = kappa.prob))
}

