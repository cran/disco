.check_matched_ortholog <- function(object) {
   # XXX
}


#' S4 class for disco
#' 
#' S4 class for disco
#' 
#' Class for holding matched orthologs and associated data. See
#' \code{makeMatchedOrholog} for details.
#' @param object An object of the matchedOrtholog class
#' @name matchedOrtholog-class
#' @rdname matchedOrtholog-class
#' @seealso makeMatchedOrholog
#' @importFrom methods setClass setMethod loadMethod is new
NULL

#' @rdname matchedOrtholog-class
setClass("matchedOrtholog", slots=list(.Data="list"), validity=.check_matched_ortholog)

# #' @rdname matchedOrtholog-class
# setMethod("[", signature(x="matchedOrtholog", i="ANY")) {
#
#   function(x, i) {
#
#     .data <- x$.Data
#     sel <- match(i, .data$
#
#
#
#
#   }
#
#
# }

#' @rdname matchedOrtholog-class
setMethod("show", "matchedOrtholog", 
  function(object) {
    n <- length(object$genes)
    cat(sprintf("Object of class matchedOrtholog\n\tdata sets: %s, %s\ngene pairs:%d\n", object$names[1], object$names[2], n))
  })

#' @name as
#' @rdname matchedOrtholog-class
setAs("matchedOrtholog", "data.frame", function(from) {
  nn <- from$names

  lfc <- from$lfc 
  colnames(lfc) <- paste0("lfc_", from$names)
  pval <- from$pval
  colnames(pval) <- paste0("pval_", from$names)

  ret <- data.frame(gene=from$genes, cbind(lfc, pval, from$extra))
  if(!is.null(rownames(from$lfc))) rownames(ret) <- rownames(from$lfc)
  return(ret)
  })
  

#' Create an object of class matchedOrtholog
#'
#' Create an object representing regulation in two sets of genes (e.g. from
#' two organisms). 
#' 
#' The object contains a mapping between two sets of genes and the
#' associated effect sizes (i.e. log fold changes) and p-values used to
#' calculate disco score and other statistics.
#' @param names character vector with names of the matched sets
#' @param genes names associated with the ortholog pairs
#' @param lfc1,lfc2 log fold changes or effect sizes
#' @param pval1,pval2 vectors of p values for the two data sets
#' @param row.names character vector to be used as row names for the data set
#' @param extra a data frame or a vector with the same length or number of rows as "genes"
#' @seealso matchedOrtholog-class
#' @export
makeMatchedOrtholog <- function(names, genes, lfc1, lfc2, pval1, pval2, row.names=NULL, extra=NULL) {
  if(!is.character(genes)) stop("genes must be a character vector")
  if(!all(c(length(lfc1), length(lfc2), length(pval1), length(pval2)) == length(genes)))
    stop( "all objects: genes, pval1, pval2, lfc1, lfc2 must have the same length")
  if(length(names) != 2) stop("names must be a character vector with length 2")

  lfc <- cbind(lfc1, lfc2)
  colnames(lfc) <- names

  pval <- cbind(pval1, pval2)
  colnames(pval) <- names

  if(!is.null(row.names)) names(genes) <- rownames(lfc) <- rownames(pval) <- row.names

  ret <- list(names=names, genes=genes, lfc=lfc, pval=pval, extra=extra)
  new("matchedOrtholog", ret)
}



