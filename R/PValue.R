###function to compute adjusted p-values

pvaluef <- function(data, mod) {
    u.gene = unique(data$Gene)
    u.gene = u.gene[order(u.gene)]
    coeftrt = pvalue = result = array(0, c(1, length(u.gene)))
    for (i in 1:length(u.gene)) {
        dat = data[data$Gene == u.gene[i], ]
        lm.dat = lm(mod, data = dat)
        coeftrt[i] = coefficients(summary(lm.dat))[2, 1]
        pvalue[i] = coefficients(summary(lm.dat))[2, 4]
    }
    result = p.adjust(pvalue, method = "fdr", n = length(pvalue))
    list(result = result, pvalue = pvalue, coef = coeftrt)
}
