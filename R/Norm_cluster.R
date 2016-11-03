
## Adaptive clustering method Excluding data with Ct values larger than 25 Cut the dendrogram tree
## from bottom up according to the height of the tree.  Returns genes and cut groups



norm.cluster = function(data, dis.method = "Euclidean", hclust.method = "complete", ct = 25, missing = 0,
    size = 10) {

    s.dat = data[, c("Sample", "Gene", "Ct")]
    u.sample = unique(s.dat$Sample)
    colnames(s.dat) = c("Sample", "Gene", "value")
    dat.cast = reshape2::acast(s.dat, Gene ~ Sample)
    dat = scr(dat.cast, threshold = ct, missing.percent.thresh = missing)
    u.gene = rownames(dat)
    imp.dat = impmat(dat)

    st.dat <- t(scale(t(imp.dat)))


    if (dis.method == "Euclidean") {
        center.dat = t(scale(t(imp.dat), scale = FALSE))
        dat.dist = dist(center.dat)
    } else if (dis.method == "1-Cor") {
        dat.dist = as.dist((1 - cor(t(st.dat)))/2)
    }

    dat.hclust = hclust(dat.dist, method = hclust.method)
    height = dat.hclust$height

    k = size + 1
    if (k < length(height)) {
        group = cutree(dat.hclust, h = height[k])
    } else {
        k = length(height)
        group = cutree(dat.hclust, h = height[k])
    }

    while (max(table(group)) < size && k < length(height)) {
        k = k + 1
        group = cutree(dat.hclust, h = height[k])
    }
    groups = which(table(group) >= size)

    ## if none of the groups are greater than size, take the largest available group
    if (length(groups) == 0) {
        groups = which(table(group) == max(table(group)))
    }

    if (length(groups) == 1)
        selected = groups else {
        dist = array(0, c(length(groups), 1))
        for (i in 1:length(groups)) {
            gene.clust.test = names(which(group == groups[i]))
            mat.clust.test = array(0, c(length(gene.clust.test), length(u.sample)))
            mat.clust.test = dat[gene.clust.test,]
            dist[i] = dist.func(mat.clust.test, standardize = FALSE)
        }
        selected = groups[which.min(dist)]
    }

    gene.clust = names(which(group == selected))



    return(list(gene.clust = as.character(gene.clust), cut.group = group))

}



