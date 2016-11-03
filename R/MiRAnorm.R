
#' Adaptive algorithm to identify normalization genes.
#'
#' \code{miranorm} returns a list of suggested normalizaiton genes based on supplied data.  Various output figures are also available
#' if requested in the parameter list.
#'
#' @param data Dataframe containing at a minimum:  Sample, Gene, Ct, and Trt
#' @param group Treatment allocation.  This should be the same length as the number of rows in data.
#' @param dis.method Distance metric used to calculate pairwise distance between individual miRNAs across all samples.
#'        Currently "Euclidean" and "1-Cor" are implemented.
#' @param hclust.method The agglomeration method used to group genes.
#'        Methods are the same as defined in the hclust function in the stats package and include
#'        "single", "average", "complete", and "ward.D2".
#'        "single" is recommended as it is more robust to small perturbations and tends to form the "chaining" phenomenon useful for defining normalizing genes.
#' @param ct Cycle threshold values at or above this level are treated as NA for the purposes of determining normalization genes.  Recommend the value be set to 25.
#' @param missing Defines maximum percentage of samples missing for a given gene before that gene is excluded from dataset during normalization.
#' @param method Choice of "simple" or "complex".
#'        Simple runs a single pass of miranorm for suggested normalizing genes.
#'        Complex runs bootstrap samples across a range of sizes to compute a stability metric.
#' @param min When method is chosen as "complex", this determines the minimum selected size at which the stability evaluation is done.
#' @param max When method is chosen as "complex", this determines the maximum selected size at which the stability evaluation is done.
#' @param selected How many adaptive normalizing genes to search for in panel.  Note, actual number of genes found may be larger based on tree cut.
#' @param ggplot "True" or "False" to output general raw data plots.
#' @param clustplot "True" or "False" to output or suppress stability plot.  Only applicable if method = "complex".
#' @param heatmap "True or "False" to output or suppress heatmap plot.
#' @param known.positives Names of miRNA that are known positive.  These will be added automatically to the heatmap plot.
#' @param suggested.list  Names of miRNA that are user suggested normalizing miRNA.  These will be added automatically to the heatmap plot.
#' @param exclude List of miRNA to exclude from the selection process for HK genes, eg: known.positives should be included here.
#' @return A list including the following: nmean, nmed, nlcv, lcv.gene, ncls, ncls.gene.
#'
#' nmean is the dataset normalized to the global mean.
#'
#' nmed is the dataset normalized to the global median.
#'
#' nlcv is the dataset normalized to the average of the 3 genes with the lowest CV.
#'
#' lcv.gene is the names of the 3 genes with used to normalize nlcv
#'
#' ncls is the dataset normalized to the adaptive normalizing genes chosen from miranorm.
#'
#' ncls.gene is the names of the genes chosen as adaptive normalizing genes chosen from miranorm.
#'
#'
#' @examples
#' dat = simData(n.trt=15, n.ctrl=15, n.gene=30, n.err=10, sigma.error = c(1, 0.3), mean.sample = 2,
#' sigma.sample = 1.88 , sigma.gene = 0.1, n.big.effect = 5, n.small.effect = 10, mean.big.effect = 2,
#' mean.small.effect = 1.2)$sim
#'
#' obj = miranorm(data = dat, group = dat$Group, method="simple")
#'
#'
#' @export
#' @import ggplot2
#' @import grDevices
#' @import stats
#' @import graphics


## Adaptive clustering method
## Excluding data with Ct values larger than 25
## Cut the dendrogram tree from bottom up according to the height of the tree
## MiRAnorm method is compared with mean, median, lowestCV, and prespecified HK (if any)genes
## trace.selection.kappa is the main function for MiRAnorm, at this moment it has not made to be full automatic. It will plot the stability curve for MiRAnorm across different cluster levels and print out the "Stability" summary statistic. But you still need to specify how many miRNAs you are going to include for HK miRNAs after you see the curve and "Stability" summary statistic. You need to specify the number in "selected = N". The default is 4. Then it will select the top number N miRNAs with the highest stability as HK miRNAs for normalization. Refer to trace.selection.kappa function for more details about parameters.

## GGplot plots the Ct values along all the samples or a random set of the samples
## each thick line indicates the normalization factor derived from different methods
## Heatmap & Dendrogram plots clustered Ct values across the set of treatment group and ctrl group
## Volcano plot is also included



## Main function
miranorm = function(data = dat, group = dat$Trt, max = 15, min = 3, method = "complex", dis.method = "Euclidean", hclust.method = "single", ct = 25, missing = 0, clustplot = TRUE, selected = 4, ggplot = TRUE, heatmap = TRUE, known.positives = NULL, suggested.list = NULL, exclude=NULL){

   ##included for R cmd check of undefined global variables, note: these are actually defined names
  ##since the datasets intended for the package are well defined and structured.
  Col = Ct = Gene = Sample = Significance = Size = Stability = effect.size = p.value = NULL;

  dat = data
  dat$Trt = group
  volcano = FALSE ##disabled feature for this version of the package.

  ##error checking##
  Check <- ArgumentCheck::newArgCheck()

  if ( !Reduce("&", c("Sample", "Gene", "Ct") %in% names(dat)) ) {
    ArgumentCheck::addError(
      msg = "'Sample', 'Gene', or 'Ct' not in data",
      argcheck = Check
    )
  }

  if ( Reduce("&", c("Sample", "Gene", "Ct")  %in% names(dat)) ) {
    ##error in Sample naming structure
    if ( sum( grepl("[.]", dat$Sample) ) > 0 ) {
      ArgumentCheck::addError(
        msg = "Unacceptable naming format for Sample.  Please refrain from using '.' in the name.",
        argcheck = Check
      )
    }

    ##error in Trt naming structure
    if ( sum( grepl("[.]", dat$Trt) ) > 0 ) {
      ArgumentCheck::addError(
        msg = "Unacceptable naming format for Trt.  Please refrain from using '.' in the name.",
        argcheck = Check
      )
    }

    ##check for unique sample/gene combinations
    tmp <- table(dat$Sample, dat$Gene) %in% c(0,1)
    if ( !all(tmp) ) {
      ArgumentCheck::addError(
        msg = "Sample*Gene combinations are not unique in the data!",
        argcheck = Check
      )
    }

    ##warn about amount of genes missing
    s.dat = data[, c("Sample", "Gene", "Ct")]
    colnames(s.dat) = c("Sample", "Gene", "value")
    dat.cast = reshape2::acast(s.dat, Gene ~ Sample)
    dat.check = scr(dat.cast, threshold = ct, missing.percent.thresh = missing)

    if (nrow(dat.check)/nrow(dat.cast) < 0.5) {
      ArgumentCheck::addWarning(
        msg = paste("Only", nrow(dat.check), "of", nrow(dat.cast), "genes meet maximum ct and missing settings."),
        argcheck = Check
      )
    }
  }

  ArgumentCheck::finishArgCheck(Check)

  ###start function###
  nraw = dat

  dat$Sample <- with(dat, paste(Trt, Sample, sep="."))

  if (method == "simple"){
  	obj <- norm.cluster(dat[!dat$Gene %in% exclude, ], dis.method = dis.method,  hclust.method = hclust.method, ct = ct, missing = missing, size = selected)
  	gene.clust <- obj$gene.clust;
  	kappa.prob=0;
  } else {
  n.core <- parallel::detectCores() - 1
  obj = trace.selection.kappa(dat[!dat$Gene %in% exclude, ], id.var="Sample", n.boot=200, max.clust=max, min.clust = min, dis.method = dis.method, hclust.method = hclust.method, plot = FALSE, selected = selected, multicore = n.core)
  co.mat <- obj$co
  kappa.prob <- obj$kappa; ##agreement##
  co.mean <- apply(co.mat, 1, mean)
  co.order <- co.mean[order(co.mean, decreasing = TRUE)]
  print(co.order[1:10])
  gene.clust <- names(co.order[1:selected])
  }
  ndat = norm.benchmark(dat, ct = ct, missing = missing, gene.clust = gene.clust)
  nmean = ndat$normal.mean
  nmed = ndat$normal.median
  nlcv= ndat$normal.hkg
  gene.clust0 = ndat$u.gene
  ncls = ndat$normal.clust

  if (clustplot == TRUE & method == "complex"){

    pdf("Clustplot.pdf",width=15,height=10)
    co.long <- reshape2::melt(co.mat, id.vars = rownames(co.mat), measure.vars = colnames(co.mat), variable.name = "Stability")
    colnames(co.long) = c("Gene", "Size", "Stability")
    co.long$Col = "Non_selected";
    co.long$Col[co.long$Gene %in% names(co.order[1:selected])] <- as.character(co.long$Gene[co.long$Gene %in% names(co.order[1:selected])])
    col.vec = c( "dark blue", "blue", "dark red", "red", "brown", "purple", "pink", "orange", "yellow" ,"dark green", "green", "gold", "black")
    print(ggplot(co.long, aes(x = Size, y = Stability, group = Gene, col = Col)) + geom_line(alpha = 0.5) + labs(title="Stability") + theme_bw()+
      scale_colour_manual(values = c(col.vec[1:selected], "grey")))
    dev.off()

  }

  if (ggplot == TRUE){
      dat74 = dat
      dat.wide <- reshape(dat74[, c("Sample", "Gene", "Ct")], timevar ="Sample", idvar= "Gene", direction="wide")
      names(dat.wide)[2:ncol(dat.wide)] <- sub("Ct[.]", "", names(dat.wide)[2:ncol(dat.wide)])
      u.gene <- unique(dat74$Gene)
      ##length(u.gene)
      u.sample = unique(dat74$Sample)

      name.gene = array(0, c(length(u.gene), 2))
      colnames(name.gene) = c("Gene.Name", "CV.Rank")
      dat.wide$cv = apply(dat.wide[, 2:ncol(dat.wide)], 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))

      dat.temp <- dat.wide[dat.wide$Gene %in% gene.clust,]
      dat.temp.lcv <- dat.wide[dat.wide$Gene %in% ndat$lcv,]


      #####
      datplot1 <- datplot2 <- datplot3 <- dat74[c("Sample", "Gene", "Ct")];
      datplot1$Method = "A. MiRAnorm"; datplot2$Method = "B. Global Mean"; datplot3$Method = "C. LowestCV";

      datplot <- rbind(datplot1, datplot2, datplot3)

      ##miranorm panel
      datplot.MiRAnormIndiv <- datplot1[datplot1$Gene %in% gene.clust,]
      datplot.MiRAnormAvg <- data.frame(Sample=names(dat.temp)[2:(ncol(dat.temp)-1)], Gene="HK.avg",
                        Ct=colMeans(dat.temp[,2:(ncol(dat.temp)-1)], na.rm = TRUE), Method = "A. MiRAnorm")

      ##Global Panel
      datplot.GlobalAvg <- data.frame(Sample=names(dat.wide)[2:(ncol(dat.wide)-1)], Gene="GlobalMean",
                                      Ct=colMeans(dat.wide[,2:(ncol(dat.wide)-1)], na.rm = TRUE), Method = "B. Global Mean")

      #Ct.median = apply(dat.wide[, 2:(ncol(dat.wide)-1)], 2, median, na.rm = TRUE)
      #datplot.GlobalMed <- data.frame(Sample=names(dat.wide)[2:(ncol(dat.wide)-1)], Gene="GlobalMedian", Ct=Ct.median, Method = "Global Mean")

      ##Lowest CV Panel
      datplot.LCVIndiv <- datplot3[datplot3$Gene %in% ndat$lcv,]
      datplot.LCVAvg <- data.frame( Sample=names(dat.temp.lcv)[2:(ncol(dat.temp.lcv)-1)], Gene="LCV.avg",
                                    Ct=colMeans(dat.temp.lcv[,2:(ncol(dat.temp.lcv)-1)], na.rm = TRUE), Method = "C. LowestCV")

      if ( length(suggested.list) != 0 ) {
        datplot4 <-  dat74[c("Sample", "Gene", "Ct")];
        datplot4$Method = "D. Endogenous Control";
        datplot <- rbind(datplot1, datplot2, datplot3, datplot4)

        ##Endogenous Control Panel
        datplot.EndoIndiv <- datplot4[datplot4$Gene %in% suggested.list,]
      }

      ##suppress ggplot warnings
      oldw <- getOption("warn")
      options(warn = -1)
      if (length(suggested.list)==0){
        pdf("GGplots.pdf",width=10,height=9)
        print(ggplot(datplot, aes(x=Sample, y=Ct, group=Gene)) + geom_line(alpha=0.5, col="grey") + facet_wrap(~Method, ncol=1) +
        geom_line(data=datplot.MiRAnormIndiv, aes(x=Sample, y=Ct, group=Gene), col="darkblue") +
        geom_line(data=datplot.MiRAnormAvg, aes(x=Sample, y=Ct, group=Gene), col="blue", lwd=2, alpha=0.3) +
        geom_line(data=datplot.GlobalAvg, aes(x=Sample, y=Ct, group=Gene), col="red", lwd=2, alpha=0.3) +
        geom_line(data=datplot.LCVIndiv, aes(x=Sample, y=Ct, group=Gene), col="darkgreen") +
        geom_line(data=datplot.LCVAvg, aes(x=Sample, y=Ct, group=Gene), col="green", lwd=2, alpha=0.3) + theme_bw() +
        theme(axis.text.x=element_text(angle = -90, hjust = 0)))
        dev.off()
      } else {
        pdf("GGplots.pdf",width=10,height=12)
        print(ggplot(datplot, aes(x=Sample, y=Ct, group=Gene)) + geom_line(alpha=0.5, col="grey") + facet_wrap(~Method, ncol=1) +
                geom_line(data=datplot.MiRAnormIndiv, aes(x=Sample, y=Ct, group=Gene), col="darkblue") +
                geom_line(data=datplot.MiRAnormAvg, aes(x=Sample, y=Ct, group=Gene), col="blue", lwd=2, alpha=0.3) +
                geom_line(data=datplot.GlobalAvg, aes(x=Sample, y=Ct, group=Gene), col="red", lwd=2, alpha=0.3) +
                geom_line(data=datplot.LCVIndiv, aes(x=Sample, y=Ct, group=Gene), col="darkgreen") +
                geom_line(data=datplot.LCVAvg, aes(x=Sample, y=Ct, group=Gene), col="green", lwd=2, alpha=0.3) +
                geom_line(data=datplot.EndoIndiv, aes(x=Sample, y=Ct, group=Gene), col="magenta") +
                theme_bw() + theme(axis.text.x=element_text(angle = -90, hjust = 0)))
        dev.off()
      }
      options(warn = oldw)
  }
    if (heatmap == T){
      dat74 = dat[dat$Gene %in% gene.clust0, ]

      dat.heat <- dat74[, c("Sample", "Gene", "Ct")]
      dat.wide <- reshape(dat.heat, timevar ="Sample", idvar= "Gene", direction="wide")


      if (dis.method == "Euclidean")  {
          toplot <- t(scale(t(dat.wide[,-1]), scale = FALSE)); ## scale the dataset by each gene across all the samples
          rownames(toplot) <- dat.wide$Gene
          dat.dist = dist(toplot)
          } else if (dis.method == "1-Cor") {
          toplot <- t(scale(t(dat.wide[,-1]))); ## scale the dataset by each gene across all the samples
          rownames(toplot) <- dat.wide$Gene
          dat.dist = as.dist((1-cor(t(toplot)))/2)
      }

      hcG.dcor <-  hclust(dat.dist, method = hclust.method)

      dat.distS = dist(t(toplot))
      hcS.dcor <- hclust(dat.distS, method = "ward.D2")

      #color pattern for Genes in columns by selected from miranorm vs not
      x.col = rep("black", length(dat.wide$Gene))
      names(x.col) = as.character(dat.wide$Gene)
      x.col[dat.wide$Gene %in% gene.clust] = "blue"

      #color pattern for Samples in rows by treatment
      u.trt <- as.character(unique(dat74$Trt))

      ##color list
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length=n+1)
        hcl(h=hues, l=65, c=100)[1:n]
      }

      col.list <- gg_color_hue(length(u.trt))

      x.col.S = rep("black", ncol(toplot))

      for (r in 1:length(u.trt)) {
        x.col.S[grep(u.trt[r], colnames(toplot))] = col.list[r]
      }

      #x.col.S[grep("Trt", colnames(toplot))] = "pink"
      #x.col.S[grep("Ctrl", colnames(toplot))] = "dark green"

      if (length(suggested.list) != 0){
        x.col1 = rep("black", length(dat.wide$Gene))
        names(x.col1) = as.character(dat.wide$Gene)
        x.col1[dat.wide$Gene %in% suggested.list] = "purple"
      }

      if (length(known.positives) != 0){
        x.col2 = rep("black", length(dat.wide$Gene))
        names(x.col2) = as.character(dat.wide$Gene)
        x.col2[dat.wide$Gene %in% known.positives] = "magenta"
      }


      if (length(known.positives)!=0 && length(suggested.list)!=0){
        pdf("Heatmap.pdf", width=10, height=10)
        heatmap3(x=data.matrix(t(toplot)),
                 Rowv = as.dendrogram(hcS.dcor) , Colv = as.dendrogram(hcG.dcor),
                 ColSideColors= rbind("AHK"=x.col, "Endo" = x.col1, "G+" = x.col2),
                 RowSideColors= rbind("Subj"= x.col.S), labCol=rownames(toplot),
                 labRow=NA,high="green", low="red",
                 scale="none",
                 na.rm = F, cexCol=0.2 + 1/log10(ncol(t(toplot))),margins = c(5, 13), zlm=c(-2,2), main = "Heatmap & Dendrogram",cex.main=1,
                 xlab = NULL, ylab = NULL)
        legend("topright",legend=c( "Known Positive",  "Prespecified Genes", "Adaptive HK Genes", u.trt), fill=c("magenta", "purple", "blue", col.list), border=FALSE, bty="n", y.intersp=0.8, cex=0.6)
        dev.off()
      } else if (length(known.positives)!=0 && length(suggested.list)==0){
          pdf("Heatmap.pdf", width=10, height=10)
          heatmap3(x=data.matrix(t(toplot)),
                 Rowv = as.dendrogram(hcS.dcor) , Colv = as.dendrogram(hcG.dcor),
                 ColSideColors= rbind("AHK"=x.col, "G+" = x.col2),
                 RowSideColors= rbind("Subj"= x.col.S), labCol=rownames(toplot),
                 labRow=NA,high="green", low="red",
                 scale="none",
                 na.rm = F, cexCol=0.2 + 1/log10(ncol(t(toplot))),margins = c(5, 13), zlm=c(-2,2), main = "Heatmap & Dendrogram",cex.main=1,
                 xlab = NULL, ylab = NULL)
          legend("topright",legend=c( "Known Positive",  "Adaptive HK Genes", u.trt), fill=c("magenta", "blue", col.list), border=FALSE, bty="n", y.intersp=0.8, cex=0.6)
          dev.off()
      } else if (length(known.positives)==0 && length(suggested.list)!=0){
        pdf("Heatmap.pdf", width=10, height=10)
        heatmap3(x=data.matrix(t(toplot)),
                 Rowv = as.dendrogram(hcS.dcor) , Colv = as.dendrogram(hcG.dcor),
                 ColSideColors= rbind("AHK"=x.col, "Endo" = x.col1),
                 RowSideColors= rbind("Subj"= x.col.S), labCol=rownames(toplot),
                 labRow=NA,high="green", low="red",
                 scale="none",
                 na.rm = F, cexCol=0.2 + 1/log10(ncol(t(toplot))),margins = c(5, 13), zlm=c(-2,2), main = "Heatmap & Dendrogram",cex.main=1,
                 xlab = NULL, ylab = NULL)
        legend("topright",legend=c( "Prespecified Genes",  "Adaptive HK Genes", u.trt), fill=c("purple", "blue", col.list), border=FALSE, bty="n", y.intersp=0.8, cex=0.6)
        dev.off()
      } else{
        pdf("Heatmap.pdf", width=10, height=10)
        heatmap3(x=data.matrix(t(toplot)),
                 Rowv = as.dendrogram(hcS.dcor), Colv = as.dendrogram(hcG.dcor),
                 ColSideColors= rbind("AHK"=x.col),
                 RowSideColors= rbind("Subj"= x.col.S), labCol=rownames(toplot),
                 labRow=NA,high="green", low="red",
                 scale="none",
                 na.rm = F, cexCol=0.2 + 1/log10(ncol(t(toplot))),margins = c(5, 13), zlm=c(-2,2), main = "Heatmap & Dendrogram",cex.main=1,
                 xlab = NULL, ylab = NULL)
        legend("topright",legend=c("Adaptive HK Genes", u.trt), fill= c("blue", col.list), border=FALSE, bty="n", y.intersp=0.8, cex=0.6)
        dev.off()
      }
    }







    if (volcano == TRUE){

      HK.dat <- ncls
      HK.dat$Gene <- factor(as.character(HK.dat$Gene))

      t.test.sum <- function(d){
        effsize <- mean(d$Ct[which(d$Trt == 1)], na.rm=T) - mean(d$Ct[which(d$Trt == 0)], na.rm=T)
        if(sum(!is.na(d$Ct[d$Trt == 0])) > 1 && sum(!is.na(d$Ct[d$Trt == 1])) > 1){
          test <- t.test(Ct ~ Trt, d)
          return(c(effsize, test$p.value))
        } else {return(c(effsize, NA))}
      }

      hk.results <- plyr::ddply(HK.dat, plyr::.(Gene), t.test.sum)
      names(hk.results)[2:3] <- c("effect.size","p.value")
      hk.results$adj.p.value <- p.adjust(hk.results$p.value)

      m.dat <- nmean
      m.dat$Gene <- factor(as.character(m.dat$Gene))

      # Calculate effect sizes and p-values for treatment v control, for each gene
      m.results <- plyr::ddply(m.dat, plyr::.(Gene), t.test.sum)
      names(m.results)[2:3] <- c("effect.size","p.value")
      m.results$adj.p.value <- p.adjust(m.results$p.value)

      # Volcano plots
      hk.results <- data.frame(hk.results, "HK")
      m.results <- data.frame(m.results, "Overall mean")
      names(hk.results)[5] <- names(m.results)[5] <- "norm.method"
      results <- rbind(hk.results, m.results)

      results$Significance <- NA
      results$Significance[which(results$adj.p.value < 0.05)] <- "significant (adj p-val < 0.05)"
      results$Significance[which(results$adj.p.value >= 0.05 & results$p.value < 0.05)] <- "significant (p-val < 0.05)"
      results$Significance[which(results$p.value >= 0.05)] <- "not significant"
      results$Significance <- factor(results$Significance) #adj.p.value < 0.05, labels=c("not significant","significant (adj p-val < 0.05)"))

      res.lab <- results[which(results$Gene %in% gene.clust),]


      p1 <- ggplot(results, aes(x=effect.size,y=-log10(p.value),col=Significance)) + geom_point()
      p2 <- geom_vline(xintercept=0, linetype="dashed")
      p3 <- geom_text(data=res.lab, aes(x=effect.size,y=-log10(p.value),label=Gene), size=2, hjust=0, show.legend=F)
      p1+p2+p3+facet_wrap(~norm.method) + labs(x="Effect size (Ct)",y="-log p-value") + theme_bw() + theme(strip.text.x=element_text(size=10),strip.text.y=element_text(size=8),plot.title=element_text(size=12),axis.title=element_text(size = 12),axis.text=element_text(size=10),legend.title=element_text(size=10),legend.text=element_text(size=10))#, axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
      ggsave("volcano.png",width=14,height=6)

    }

  ##check gene.clust and issue warning if there are apparent treatment effects
  if (nlevels(nraw$Trt) >= 2) {

    dat.analy <- nraw[nraw$Gene %in% gene.clust, ]
    dat.analy$Gene <- factor(dat.analy$Gene)

    trt.mod <- function(d){
      mod <- aov(data=d, Ct ~ Trt)
      summary(mod)[[1]][["Pr(>F)"]][[1]]
    }

    res.table = plyr::ddply(dat.analy, plyr::.(Gene), trt.mod)
    res.table$p.adj = p.adjust(res.table[,2], method="BH") #adjusted p.values

    if ( any(res.table$p.adj < 0.05 )   ) {
      msg = paste(paste(res.table[,1][res.table[,3] < 0.05  ], collapse=", "), "may have significant treatment effect.")
      warning(msg)
    }

    #k=1
    #dat.analy <- nraw[nraw$Gene == gene.clust[k],]
    #mod <- aov(data=dat.analy, Ct ~ Trt)
    #summary(mod)[[1]][["Pr(>F)"]][[1]]


  }


      # if (length(known.positives)!=0){
      #   trtset = known.positives
      #   ctrlset = u.gene[-which(u.gene %in% trtset)]
      #   dat.list = vector("list", length = 5)
      #   sim_result = vector("list", length = 5)
      #   effect = effect1 = name = array(0, c(length(u.gene), 5))
      #   mod = "Ct ~ Group"
      #   dat.list[[1]] = nraw
      #   dat.list[[2]] = nmean
      #   dat.list[[3]] = nmed
      #   dat.list[[4]] = nlcv
      #   dat.list[[5]] = ncls
      #   for (i in 1:5){
      #     obj = pvaluef(dat.list[[i]], mod)
      #     es = obj$coef
      #     effect[, i] = ifelse(abs(es)>1, 1, 0)
      #     fdr = obj$result
      #     name[, i] = ifelse(fdr<=0.1, 1, 0)
      #     sim_result[[i]] = array(0, c(4, 1))
      #     sim_result[[i]][1] = sum(name[u.gene %in% trtset, i]!=0)
      #     sim_result[[i]][2] = sum((effect[u.gene %in% trtset, i]+name[u.gene %in% trtset, i])==2)
      #     sim_result[[i]][3] = sum(name[u.gene %in% ctrlset, i]!=0)
      #     sim_result[[i]][4] = sum((effect[u.gene %in% ctrlset, i]+name[u.gene %in% ctrlset, i])==2)
      #   }
      # return(list(nmean = nmean, nmed = nmed, nlcv = nlcv, lcv.gene = ndat$lcv, ncls = ncls, ncls.gene = gene.clust, result = sim_result))
      # }

  return(list(nmean = nmean, nmed = nmed, nlcv = nlcv, lcv.gene = ndat$lcv, ncls = ncls, ncls.gene = gene.clust, kappa.prob=kappa.prob))
}


