times_to_elen <- function(phy, ts)
{
    n <- length(phy$edge.length)
    phy$edge.length <- sapply(1:n, function(i) ts[phy$edge[i,1]] - ts[phy$edge[i,2]])
    return(phy)
}

sample2tree <-function(x, s_idx)
{
    m <- length(x$undated$tip.label)
    n <- 2*m-1

    tiplabs <- x$undated$tip.label
    edge <- matrix(-1L, n, 2)

    time_cols <- paste0("t_", 1:n)
    pa_cols <- paste0("pa_", 1:n)

    ts <- suppressWarnings(unlist(x$draws[s_idx, time_cols]))
    pas <- suppressWarnings(unlist(x$draws[s_idx, pa_cols]))

    edge[,2] <- 1L:n
    edge[,1] <- pas
    edge <- edge[-which(is.na(edge[,1])),]

    elen <- sapply(1:(n-1), function(i) abs(ts[edge[i,1]] - ts[edge[i,2]]))

    tr <- list(edge=edge, edge.length=elen, Nnode=m-1L, tip.label=tiplabs)
    class(tr) <- "phylo"

    return(tr)
}

maj_rule_tree <- function(trees)
{
    tip_labs <- trees[[1]]$tip.label
    edge <- maj_rule_edge(trees)

    k <- nrow(edge)
    m <- length(tip_labs)
    elen <- rep(1, k)

    tr <- list(edge=edge, edge.length=elen, Nnode=k+1-m, tip.label=tip_labs)
    class(tr) <- "phylo"
    return(tr)
}
