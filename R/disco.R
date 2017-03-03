.pkgenv <- new.env(parent=emptyenv())

.getgenes <- function(x, g=NULL) {
  if(is.null(g)) return(list(all=x$genes))

  if(is.vector(g)) {
    g <- list(unnamed=g)
  } else if(class(g) == "tmod") {
    g <- g$MODULES2GENES
  } else {
    stop("unsupported class for object g")
  }
  g
}






## ----fun 1---------------------------------------------------------------

#' Calculate concordance score for differentially expressed genes in two data sets
#' 
#' For each pair of matched (orthologous) genes from an object of class mtchedOrtholog calculates the concordance / discordance score.
#' The rationale is to formulate a metric allowing to compare heterologous
#' data sets, e.g. data sets from different sources, technical platforms or
#' even organisms.
#' 
#' \deqn{\textrm{disco.score} := |logFC_1|\cdot|logFC_2|\cdot(|log_{10}(p_1)| + |log_{10}(p_2)|) \cdot \textrm{sign}(logFC_1) \cdot \textrm{sign}(logFC_2)}
#'
#' Where
#' \describe{
#'    \item{\eqn{logFC_1}}{log fold change of the expression change of the gene from the first data set}
#'    \item{\eqn{logFC_2}}{log fold change of the expression change of the gene from the second data set}
#'    \item{\eqn{p1}}{differential regulation p-value for the gene from the first data set}
#'    \item{\eqn{p2}}{differential regulation p-value for the gene from the second data set}
#' }
#'
#' @param x object of class matchedOrtholog
#' @return A numerical vector of the disco scores
#' @examples 
#' data(orthologs)
#' orthologs$disco.score <- disco.score(orthologs)
#' @export
disco.score <- function (x) {
  logFC1 <- x$lfc[,1]
  logFC2 <- x$lfc[,2]
  p1 <- x$pval[,1]
  p2 <- x$pval[,2]
  if (!all(c(length(logFC1), length(logFC2), length(p1)) == length(p2)))
  stop("all vectors must have the same length")
  ret <- abs(logFC1) * abs(logFC2) * (abs(log10(p1)) + abs(log10(p2))) * sign(logFC1) * sign(logFC2)
  return(ret)
}

## ----fun 2---------------------------------------------------------------

#' Calculate correlation between regulation of orthologous genes in a chosen gene set
#' 
#' Calculates correlation between regulation of group of genes in two data sets. 
#' Takes as arguments an object of class matchedOrtholog, a vector of genes of interest
#' (for example, genes belonging to an expression module) and method of correlation calculation. The
#' function also checks how many of the genes in the vector are present in the matchedOrtholog object.
#' frame.
#'  
#' @param x an object of class matchedOrtholog
#' @param g either a vector of gene names of genes of interest, 
#'          or an object of class "tmod" with definitions of modules. If
#'          NULL, all genes will be used
#' @param method method for calculation of correlation (default: pearson)
#' @return data frame with two columns: Cor - calculated correlation between the genes, 
#'         P.Value - p value from test of significance of the correlation,
#'         NGenes - number of genes included in correlation calculation. 
#'         Row names of the data frame correspond to the names of the modules.
#' @examples 
#' library(tmod)
#' data(tmod)
#' data(orthologs)
#'
#' genes <- tmod$MODULES2GENES[["LI.M0"]]
#' a <- modCor(orthologs, genes)
#'
#' # Using tmod objects directly
#' a <- modCor(orthologs, tmod[c("LI.M0", "LI.M1.0")], "spearman")
#' @export
modCor <- function (x, g=NULL, method="pearson") {
  
  #if(!is.data.frame(df)) stop("df must be a data frame")
  
  x1 <- x$lfc[,1]
  x2 <- x$lfc[,2]
  gn <- x$genes

  # figure out the genes
  g <- .getgenes(x, g)


  ret <- sapply(g, function(gset) {
    sel <- gn %in% gset
    cc <- cor.test(x1[sel], x2[sel], method=method)
    c(Cor=cc$estimate, P.Value=cc$p.value, NGenes=sum(sel))
  })

  ret <- data.frame(t(ret))

  lowN <- sum(ret$NGenes < 3)
  if(lowN > 0) warning(sprintf("less than 3 genes in %d gene sets", lowN))

  return(ret)
}

## ----fun 3---------------------------------------------------------------

#' plots log2FC of the corresponding genes against each other
#' 
#' Creates a plot of log2FC of the chosen genes in one data set against log2FC of corresponding genes in the other data set. 
#' Each dot represents a gene and is colored according to disco.score calculated for this pair of heterologous genes: 
#' the stronger is red color the more concordantly regulated is the gene pair, and the stronger is the blue color, 
#' the more discordantly regulated is the gene pair.
#' 
#' @param x an object of class matchedOrtholog 
#' @param disco.score a numerical vector of values of disco.score assigned to each orthologous gene pair, ordered according to the gene order in matchedOrtholog object
#' @param g either a vector of gene names of genes of interest
#' @return plot
#' @examples 
#' library(tmod)
#' data(tmod)
#' data(orthologs)
#' ds <- disco.score(orthologs)
#' plotDisco(orthologs, ds, tmod$MODULES2GENES[["LI.M1.0"]])
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @import tmod
#' @importFrom RColorBrewer brewer.pal
#' @export
plotDisco <- function (x, disco.score, g=NULL) {
  
  colours <- brewer.pal(7, "Spectral")
  tt <- data.frame(genes=x$genes, x$lfc, disco.score)

  if(!is.null(g)) tt <- tt [which(toupper(x$genes) %in% toupper(g)), ]
  ggplot(tt, aes(x = tt[,2], y = tt[,3], col = tt[,4] )) + 
  geom_point() +
  theme_bw(30) +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  scale_colour_gradient2(low="blue", high="red", mid="grey", space="Lab") +
  labs(x = "lfc1") +
  labs(y = "lfc2")
}

## ----fun 4---------------------------------------------------------------


#' Calculate ortholog gene correlation coefficients
#'
#' Calculate correlation between sets of matched orthologous genes using
#' either method from Takao et al. (2015) or Seok et al. (2013)
#' 
#' This function calculates Pearson Correlation, squared Pearson correlation and Spearman correlation 
#' for mouse and human logFC values for genes with significant p-values for differential regulation. 
#' As in PNAS 2013 110 (9) 3507-3512, for the calculation of squared Pearson correlation coefficient 
#' all the genes that are significantly regulated in at least one specie are taken into account, 
#' while for the calculation of Spearman's correlation coefficient only the genes significantly regulated in both species
#' as in PNAS 2015 112 (4) 1167-1172.
#'
#' The parameters can be calculated either for all genes in the
#' \code{orthologs} object, for an explicit list of genes as a character
#' vector parameter \code{g}, or for a set of \code{tmod} modules, if the
#' parameter \code{g} is an object of class \code{tmod}. 
#'
#' @section Bibliography:
#' Seok, Junhee, et al. "Genomic responses in mouse models poorly mimic human
#' inflammatory diseases." Proceedings of the National Academy of Sciences
#' 110.9 (2013): 3507-3512.
#' 
#' Takao, Keizo, and Tsuyoshi Miyakawa. "Genomic responses in mouse models
#' greatly mimic human inflammatory diseases." Proceedings of the National
#' Academy of Sciences 112.4 (2015): 1167-1172.
#' 
#' @param x an object of class matchedOrtholog 
#' @param pval p value threshold
#' @param lfc log_2 fold change threshold
#' @return corOrt returns a data frame with the following columns:
#' \describe{
#'   \item{r}{Pearson correlation coefficient r}
#'   \item{Seok}{correlation measure described by Seok et al. -- squared
#'   Pearson r correlation coefficient for genes significantly up- or
#'   downregulated in one of the two orthologous sets}
#'   \item{Seok.N}{number of gene pairs that were used to calculate the Seok
#'   coefficient}
#'   \item{rho}{Spearman correlation coefficient rho}
#'   \item{Takao}{correlation measure described by Takao et al. -- Spearman
#'   correlation coefficient for genes significantly up- or downregulated in
#'   both of the two gene orthologous sets}
#'   \item{Takao.N}{number of gene pairs that were used to calculate the Takao
#'   coefficient; smaller or equal to Seok.N}
#'   \item{N}{total number of gene pairs in the given gene set}
#' }
#' @inheritParams modCor
#' @importFrom stats cor cor.test
#' @examples 
#' library(tmod)
#' data(tmod)
#' data(orthologs)
#' a <- corOrt(orthologs)
#' disco <- disco.score(orthologs)
#' ord <- order(disco, decreasing = TRUE)
#' concordant <- tmodCERNOtest(toupper(orthologs$genes)[ord])
#' corOrt(orthologs, g=tmod[concordant$ID])
#' @seealso modCor
#' @export
corOrt <- function (x, g=NULL, pval=0.05, lfc=0) {
  tt <- data.frame(genes=x$genes, x$lfc, x$pval)
  colnames (tt) <- c("genes", "lfc1", "lfc2", "pval1", "pval2")

  # process the g parameter
  g <- .getgenes(x, g)

  ret <- sapply(g, function(gg) {
    sel <- tt$genes %in% gg
    .corOrtCalc(tt[sel,], pval, lfc)
  })
 
  ret <- data.frame(t(ret))
  return(ret)
}

## calculate the seok / takao method for a selection of genes
.corOrtCalc <- function(tt, pval, lfc) {

  lfc1 <- tt$lfc1
  lfc2 <- tt$lfc2
  p1 <- tt$pval1
  p2 <- tt$pval2
  
  Pcorr <- cor(lfc1, lfc2, method="pearson")
  
  sel1 <- (p1 < pval & abs(lfc1) > lfc) | (p2 < pval & abs(lfc2) > lfc)
  Seok.N <- sum(sel1)
  Seok <- cor(lfc1[sel1], lfc2[sel1], method="pearson")
  SeokPSquared <- Seok^2
  n <- length(sel1[which(sel1==TRUE)])
  
  Scorr <- cor(lfc1, lfc2, method="spearman")
  sel2 <- (p1 < pval & abs(lfc1) > lfc) & (p2 < pval & abs(lfc2) > lfc)
  Takao.N <- sum(sel2)
  Takao <- cor(lfc1[sel2], lfc2[sel2], method="spearman")
  n1 <- length(sel2[which(sel2==TRUE)])
  
  return (c(r=Pcorr, Seok=Seok, Seok.N=Seok.N, rho=Scorr, Takao=Takao, Takao.N=Takao.N, N=nrow(tt)))
}

## ----fun 5---------------------------------------------------------------

#' Common modules in specified comparisons  
#' 
#' Summarize numbers of enriched modules in an unrestricted number of comparisons and list all the modules that belong to intersection of specified comparisons.
#' 
#' @param ... any number of data frames created as a result of tmodCERNOtest functions
#' @return list of two elements: RowN - list, each element of it contains number of rows in every subsequent comparison, 
#' int character vector with IDs of all modules belonging to the intersection of comparisons provided as the input
#' @examples 
#' library(tmod)
#' data(tmod)
#' data(orthologs)
#' disco <- disco.score(orthologs)
#' g <- toupper(orthologs$genes)
#' concordant <- tmodCERNOtest(g[order(disco, decreasing = TRUE)])
#' discordant <- tmodCERNOtest(g[order(disco)])
#' sum1 <- makeSummary(concordant, discordant)
#' @export
makeSummary <- function (...) {

  # use the names of the variables provided for the return values
  dots <- substitute(list(...))[-1]
  d <- c(sapply(dots, deparse))
  a <- list (...)
  names(a) <- d

  b <- sapply(a, function(a1) nrow(a1), simplify=F)
  aids <- lapply(a, function(a) a$ID)
  int <- Reduce(intersect, aids)
  list(RowN=b, int=int)
}

## ----fun 6---------------------------------------------------------------

#' Calculates percentage of present genes from a specified superset of expression modules
#' 
#' This function counts the number of genes in modules specified in data.frame created by tmodCERNOtest function 
#' and returns percentage of those genes which are present in data.frame on the basis of which the enrichment was calculated, 
#' as well as original number of genes present in the module.
#' 
#' @param orthologs object of class matchedOrtholog
#' @param mset a tmod module set (object of class tmod). If NULL, the default tmod object will be used.
#' @return data frame with columns specifying percentages of genes present in subsequent modules "p_genes" 
#' and the overall amount of genes in the modules "all_genes"
#' @examples 
#' library(tmod)
#' data(tmod)
#' data(orthologs)
#' disco <- disco.score(orthologs)
#' ord <- order(disco, decreasing = TRUE)
#' concordant <- tmodCERNOtest(toupper(orthologs$genes)[ord])
#' modtable <- percentGenes(orthologs, tmod[concordant$ID])
#' cbind(concordant, modtable)
#' @importFrom utils data
#' @export
percentGenes <- function (orthologs, mset=NULL) {

  if(is.null(mset)) {
    data("tmod", package="tmod", envir=.pkgenv)
    mset <- .pkgenv$tmod
  }

  ids <- mset$MODULES$ID
  gn <- orthologs$genes
  p_genes <- sapply (ids, function(id) sum(mset$MODULES2GENES[[id]] %in% gn) / length(mset$MODULES2GENES[[id]]))
  p_genes <- as.data.frame(p_genes)
  all_genes <- sapply(ids, function(id) length(mset$MODULES2GENES[[id]]))
  enr.result <- cbind(p_genes, all_genes)
  return(enr.result)
}
