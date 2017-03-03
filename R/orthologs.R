#' Example ortholog data set
#'
#' Example object of class matchedOrtholog containing mouse and human data sets
#'
#' The matchedOrtholog object \code{orthologs} contains matched orthologous
#' genes from mouse and human samples, and the results of differential
#' expression analysis as \eqn{log_2} fold changes and p-values.  
#'
#' The human data set contains
#' blood transcriptional profiles of 46 TB patients and 
#' 62 healthy individuals and is available on the Gene Expression Omnibus
#' database under the accession number GSE28623 (Maertzdorf et al., 2011). 
#'
#' The mouse data set was derived from 129S2 mice infected with TB for 24h (5 mice) and 
#' before infection (5 mice) and is available under the accession number
#' GSE89392.  The data sets have been analyzed with limma R package for differential
#' expression analysis (Ritchie et al., 2015). 
#' The data sets were background corrected using the normexp method and
#' quantile normalized between arrays. 
#' Limma lmFit function was used to fit linear models which included the
#' factors: stimulus type and time point. 
#' The p-values were calculated based on the moderated t-statistic and most
#' differentially regulated genes were retrieved with topTable function.
#' Orthologous genes were assigned to each other between corresponding
#' human and mouse data sets used in each comparison. 
#' Probe names specific to the microarray used were assigned an ENSEMBL
#' identifier with use of mapIds function from the
#' biomaRt package (version 2.24.1, Durinck et al., 2005; Durinck,
#' Spellman, Birney, & Huber, 2009). 
#'
#' Multiple repeating probes were averaged by applying limma
#' \code{avereps} function. Then, orthologous human and mouse genes were
#' identified 
#' with biomaRt \code{getLDS} function based on homology mapping between
#' different species interlinked in Ensembl data base 
#' (with attributes and filters defined as \code{ensembl_gene_id}.  Only
#' the putative orthologs with a 1:1 mapping (no potential in-paralogs) 
#' were included in the further analysis.
#' 
#' The gene names in the object correspond to the human gene names.
#' @name orthologs
#' @examples
#' data(orthologs)
#' # view the object as a data frame
#' head(as(orthologs, "data.frame"))
#' # calculate the disco score
#' ds <- disco.score(orthologs)
NULL

