

###########################
###heatmap 3. modification of heatmap2 to allow multiple sidebars###
########################

heatmap3 <- function (x, imp = TRUE, Rowv = NA, Colv = NULL, distfun = dist,
	hclustfun = hclust, add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
	scale = c("none","row","column"), na.rm = TRUE, margins = c(5,5), ColSideColors, RowSideColors,
	cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, labCol.line=-0.5,
	totalC=nc, main = NULL, xlab = NULL, ylab = NULL, verbose = getOption("verbose"),
                      cex.main=1.5,
	methodR = "ward", methodC = "ward", zlm = c(-0.5, 0.5), high="red", low="green", mid="black",
	addamps=NULL, colamps=NULL, cexAmp=.25, ...)
{
	## 4. Ritu
	if (!is.matrix(x)) {
		x <- matrix(x,ncol=1)
		if (!is.null(addamps)) {
			addamps <- matrix(addamps,ncol=1)
		}
	}
    scale <- if (symm && missing(scale))
	"none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
	## 4. Ritu
    #if (nr <= 1 || nc <= 1) stop("`x' must have at least 2 rows and 2 columns")
    if (nr <= 1 || nc < 1) stop("`x' must have at least 2 rows and 1 column")
    if (!is.numeric(margins) || length(margins) != 2) stop("`margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x), method = methodR)
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv)
                ddr <- reorder(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc)
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm)
                x
            else t(x)), method = methodC)
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv)
                ddc <- reorder(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong\nlength")
    }
    else colInd <- 1:nc
        x <- x[rowInd, colInd]
        if (is.null(labRow)) {
          labRow <- if (is.null(rownames(x)))
            (1:nr)#yyfix remove [rowInd]
          else rownames(x)}
        else
          labRow <- labRow[rowInd]
        if (is.null(labCol))
          labCol <- if (is.null(colnames(x)))
            (1:nc)#yyfix
          else colnames(x)
        else
          labCol <- labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    #lwid <- c(if (doRdend) 1 else 0.05, 4)
    #lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 4)

	## 7. Ritu
	#lwid <- c(1, 4)
	#lhei <- c(1 + if (!is.null(main)) 0.2 else 0, 4)
	lwid <- c(ifelse(doRdend,0.5,0.01), 4)
	lhei <- c(ifelse(doCdend,0.5,0.01) + ifelse(is.null(main),0,0.2), 4)

	## 2. Ritu
	if (!is.null(ColSideColors)) {
	#if (!missing(ColSideColors)) {
       # if (!is.character(ColSideColors) || length(ColSideColors) !=
        #    nc)
          #  stop("'ColSideColors' must be a character vector of\nlength ncol(x)")

		## 6. Ritu
		if (!is.matrix(ColSideColors)) {
			ColSideColors <- matrix(ColSideColors,nrow=1)
			rownames(ColSideColors) <- ""
		}

	## 11. Ritu
	#lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
	lmat <- rbind(lmat[1, ] + 1, c(rep(NA, ncol(lmat) - 1),1), lmat[2, ] + 1)

        lhei <- c(lhei[1], 0.2, lhei[2])
    }
	## 3. Ritu
    if (!is.null(RowSideColors)) {
    #if (!missing(RowSideColors)) {
        #if (!is.character(RowSideColors) || length(RowSideColors) != nr)
        #    stop("'RowSideColors' must be a character vector of\nlength nrow(x)")
        #lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])

	## 11. Ritu
	if (!is.matrix(RowSideColors)) {
		RowSideColors <- matrix(RowSideColors,nrow=1)
		rownames(RowSideColors) <- ""
	}
	lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
	lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, "; lmat=\n")
        print(lmat)
    }

################################
#redo lmat:

	## 1. Ritu
	#if (!is.matrix(ColSideColors)) {
		#ColSideColors <- matrix(ColSideColors,nrow=1,dimnames=list(names(ColSideColors),1:nc))
	#}

	## 2. Ritu
	if (!is.null(ColSideColors)) {

	nr.my=nrow(ColSideColors)+2

	## 11. Ritu
	#nc.my=3
	#nfig=nr.my+2
	#lmat=matrix(0, nrow=nr.my, ncol=nc.my)
	#lmat[nr.my,]=c(nfig-1,1,nfig-2)
	#lmat[1, nc.my]=nfig
	#lmat[2:(nr.my-1),nc.my]=(nfig-3):2
	nc.my=2+ifelse(is.null(RowSideColors),0,nrow(RowSideColors))
	nfig=nr.my+nc.my-1
	lmat=matrix(0, nrow=nr.my, ncol=nc.my)
	if (is.null(RowSideColors)) {
		lmat[nr.my,]=c(nfig-1,nfig-2)
	} else {
		lmat[nr.my,]=c(nfig-1,1:nrow(RowSideColors),nfig-2)
		lwid <- c(ifelse(doRdend,1,0.01), rep(.2, nrow(RowSideColors)), 4)
	}
	lmat[1, nc.my]=nfig
	lmat[2:(nr.my-1),nc.my]=(nfig-3):(nc.my-1)

	## 11. Ritu
	### 10. Ritu
	#if (is.null(RowSideColors)) {
	#	lmat=lmat[,c(1,3:ncol(lmat))]
	#	lmat=lmat-1
	#	lmat[lmat==(-1)]=0
	#}

	## 7. Ritu
	#lhei=c(1, rep(.2, nrow(ColSideColors)),4)
	lhei <- c(ifelse(doCdend,0.3,0.01) + ifelse(is.null(main),0,0.2), rep(.1, nrow(ColSideColors)), 4) ### ***

	}

	notRequired=function(x) {
	## 11. Ritu
	if (!is.null(RowSideColors)) {
		nr.my=nrow(RowSideColors)+2
		nc.my=nrow(ColSideColors)+2
		nfig=nr.my+2
		lmat=matrix(0, nrow=nr.my, ncol=nc.my)
		lmat[nr.my,]=c(nfig-1,1,nfig-2)
		lmat[1, nc.my]=nfig
		lmat[2:(nr.my-1),nc.my]=(nfig-3):2
		if (is.null(ColSideColors)) {
			lmat=lmat[,c(1,3:ncol(lmat))]
			lmat=lmat-1
			lmat[lmat==(-1)]=0
		}
		lhei <- c(ifelse(doCdend,1,0.01) + ifelse(is.null(main),0,0.2), rep(.2, nrow(RowSideColors)), 4)
	}
	}


###################################
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
#cat("--- lmat widths, heights \n")
#print(lmat)
#print(lwid)
#print(lhei)
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)

#layout.show(n=nfig)
#}
#notRequired=function() {
#print(rownames(RowSideColors))
#print(rownames(ColSideColors))

	## 3. Ritu
	if (!is.null(RowSideColors)) {
    #if (!missing(RowSideColors)) {
		par(mar = c(margins[1], 0, 0, 0.5))

		## 11. Ritu
		#if (revC) {
		#	image(rbind(1:nr), col = RowSideColors[rev(rowInd)], axes = FALSE)
		#} else {
        	#	image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
		#}
		if (revC) {
			j=rev(rowInd)
		} else {
			j=rowInd
		}
		for (i in 1:nrow(RowSideColors)) {
        		image(rbind(1:nr), col = RowSideColors[i,j], axes = FALSE)
        		#image(rbind(1:nr), z=rbind(1:nr),col = RowSideColors[i,j], axes = FALSE)
        		#image(cbind(1:nr), z=cbind(1:nr),col = RowSideColors[i,j], axes = FALSE)
			mtext(side=1, text=as.character(rownames(RowSideColors)[i]), las=3, cex=1)
		}
	}
	## 2. Ritu
	if (!is.null(ColSideColors)) {
	#if (!missing(ColSideColors)) {
        par(mar = c(0.2, 0, 0, margins[2]))
	for (i in 1:nrow(ColSideColors)) {
		## 9. Ritu
		#image(cbind(1:nc), col = ColSideColors[i,colInd], axes = FALSE)
		image(cbind(1:nc), z=cbind(1:nc),col = ColSideColors[i,colInd], axes = FALSE, xlim = 0.5 + c(0, totalC))
		mtext(side=2, text=as.character(rownames(ColSideColors)[i]), las=1, cex=1)
	}
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none")
        x <- t(x)
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr
    x.floor <- x
    for (i in 1:ncol(x)) {
        ind1 <- (1:length(x[, i]))[x[, i] >= zlm[2] & !is.na(x[,i])]
        ind2 <- (1:length(x[, i]))[x[, i] <= zlm[1] & !is.na(x[,i])]
        x.floor[, i][ind1] <- rep((zlm[2] - 0.01), length(ind1))
        x.floor[, i][ind2] <- rep((zlm[1] + 0.01), length(ind2))
    }
	## 5. Ritu
    #image(1:nc, 1:nr, x.floor, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = maPalette(high = high, low = low, mid = mid), zlim = zlm, ...)
	if (nc==1) {
	    image(1:(2*nc), 1:nr, rbind(x.floor,x.floor), xlim = 0.5 + c(0, 2*totalC), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = colorRampPalette(c(low,mid,high))(256) , zlim = zlm, ...)
	} else {
	    image(1:nc, 1:nr, x.floor, xlim = 0.5 + c(0, totalC), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = colorRampPalette(c(low,mid,high))(256), zlim = zlm, ...)
	}
##################

    if (!is.null(addamps)) {
	addamps=addamps[rowInd, colInd]
	## 4. Ritu
	if (!is.matrix(addamps)) {
		addamps <- matrix(addamps,ncol=1)
	}
	for (i in 1:ncol(addamps)) {
		amp=which(addamps[,i]>0)
		## 12. Ritu
		### 8. Ritu
		##points(rep(i, length(amp)), amp, col=colamps, cex=.75, pch=20)
		#points(rep(i, length(amp)), amp, col=colamps, cex=.25, pch=20)
		points(rep(i, length(amp)), amp, col=colamps, cex=cexAmp, pch=20)
	}

     }
################
        thislas <- ifelse(nc<10,1,2)
    axis(1, 1:nc, labels = labCol, las = thislas, line = labCol.line, tick = 0, cex.axis = cexCol) #yyfix, remove [colInd]
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
	## 9. Ritu
	if (totalC==nc) {
		par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
	} else {
		par(mar = c(0, 0, if (!is.null(main)) 1 else 0, 2*(totalC-nc)+1+margins[2]))
	}
	if (doCdend) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	} else if (!is.null(main)) {
        frame()
	}
    if (!is.null(main))
        title(main, cex.main = cex.main)
    invisible(list(rowInd = rowInd, colInd = colInd))
}

