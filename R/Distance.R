#' Calculating the distance between each row of the matrix and the centroid
#'
#' \code{dist.func} returns a sum total of all distances for each row of the matrix from the centroid.
#'
#' @param mndata matrix dataset with numeric columns.
#' @param standardize is a booleon variable to denote whether columns should first be scaled to 0-1.
#' @param method dictates how distance is computed.  Currently, only "euclidean" is implemented.
#' @param centroid allows user to specify the centroid values for each column.  Otherwise, the mean of each column is taken by default.
# Calculating the euclidean distance between each row of the matrix (mndata) and the centroid

dist.func = function(mndata, standardize = TRUE, method = "euclidean", centroid = NULL) {
    # standardize == TRUE is default, with euclidean distance calculated after standardization
    if (standardize == TRUE) {
        standard.mndata = scale(mndata, center = TRUE, scale = TRUE)
    } else {
        standard.mndata = mndata
    }
    # length(centroid) == 0 is default, with centroid as the means of the rows of the mndata,
    # otherwise it can be an input data
    if (length(centroid) == 0)
        centroid = colMeans(standard.mndata, na.rm = TRUE) else centroid = centroid
    dists = array(0, c(nrow(mndata), 1))
    for (i in 1:nrow(standard.mndata)) {
        # method == 'euclidea' is default, the distance is calculated by euclidean distance
        if (method == "euclidean")
            dists[i] = dist(rbind(t(standard.mndata[i, ]), t(centroid)), method = "euclidean") else dists[i] = dist(rbind(t(standard.mndata[i, ]), t(centroid)), method = method)
    }
    # if (method == 'mahalanobis') dists = mahalanobis(standard.mndata, colMeans(standard.mndata),
    # cov.mndata)
    return(sum(dists))
}
