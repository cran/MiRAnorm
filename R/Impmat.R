#' Imputation function
#'
#'\code{impmat} imputes a dataset entered as rows=Genes, columns=Samples.   NA values are imputed by row.  Currently, impvalue function imputes these as the max value of the row, 
#' the largest observed value for a given Gene across all Samples.
#'
#' @param data dataset to be imputed.  Should contain all numeric values with rows=Genes and columns=Samples.
#' 
#' @export


## Imputation (with max, because NA's are likely to be those that have larger CT values)
impmat = function(data){

    dat.new = impvalue(1:nrow(data), data)
	return(dat.new)

}

impvalue <- function(x, currentdat){

        currentdat[x, ][which(is.na(currentdat[x, ]))] = max(currentdat[x, ], na.rm = T)
    return(currentdat)
}
