#' @export
clade_node_counts <- function(x, ...) UseMethod("clade_node_counts")

#' @export
find_clade_index <- function(x, ...) UseMethod("find_clade_index")

#' Count the number of nodes in clades
#' @param x an `ape::phylo` object
#' @param clades a list of vectors of tip labels corresponding to clades
#' @return a vector of node counts for each clade in clades
#' @export 
clade_node_counts.phylo <- function(x, clades)
{
    tiplabs <- x$tip.label
    clade_mat <- as_clade_mat(clades, tiplabs)
    nodes_in_clades(x, clade_mat)
}

#' Count the number of nodes in clades
#' @param x an `ape::multiPhylo` object
#' @param clades a list of vectors of tip labels corresponding to clades
#' @return a matrix of node counts. Rows correspond to phylogenies. Columns correspond to clades.
#' @export 
clade_node_counts.multiPhylo <- function(x, clades)
{
    tiplabs <- x[[1]]$tip.label
    clade_mat <- as_clade_mat(clades, tiplabs)
    do.call(rbind, lapply(x, function(y) nodes_in_clades(y, clade_mat)))
}

#' Find the indices of the MRCAs of clades in a phylogeny
#' @param x an `ape::phylo` object
#' @param clades a list of vectors of tip labels corresponding to clades
#' @return a vector of node indices correspoding to clade MRCAs
#' @export 
find_clade_index.phylo <- function(x, clades)
{
    tiplabs <- x$tip.label
    clade_mat <- as_clade_mat(clades, tiplabs)
    return(as.vector(find_clade_indices(list(x), clade_mat)))
}

#' Find the indices of the MRCAs of clades in a phylogeny
#' @param x an `ape::multiPhylo` object
#' @param clades a list of vectors of tip labels corresponding to clades
#' @return a matrix of clade MRCA indices. Rows correspond to phylogenies. Columns correspond to clades.
#' @export 
find_clade_index.multiPhylo <- function(x, clades)
{
    tiplabs <- x[[1]]$tip.label
    clade_mat <- as_clade_mat(clades, tiplabs)
    return(find_clade_indices(x, clade_mat))
}

as_clade_mat <- function(clades, tiplabs)
{
    return(as.matrix(do.call(rbind, lapply(clades, function(x) clade_as_vec(x, tiplabs)))))
}

clade_as_vec <- function(clade_labs, tip_labs)
{
    out <- rep(0L,length(tip_labs))
    out[tip_labs %in% clade_labs] <- 1L
    return(out)
}