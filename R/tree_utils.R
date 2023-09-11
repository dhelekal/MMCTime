coal2ape <-  function(samp_t, n_samp, coal_t, n_merge, node_name_prefix="N")
{
    x <- coal2newick(samp_t, n_samp, coal_t, n_merge, node_name_prefix) 
    return(list(phy=read.tree(text=x$str), dates=x$dates))
}

reindex_node <- function(phy, idx1, idx2) {

    phy_new <- phy
    
    ntip <- length(phy_new$tip.label)
    stopifnot("Cannot relabel a tip node"= idx1 > ntip)
    stopifnot("Cannot relabel a tip node"= idx2 > ntip )

    edge_tmp <- phy_new$edge
    if(!is.null(phy_new$node.label))
    {
        nodelabs_tmp <- phy_new$node.label
        nodelabs_tmp[idx1] <- phy_new$node.label[idx2]
        nodelabs_tmp[idx2] <- phy_new$node.label[idx1]
        phy_new$node.label <- nodelabs_tmp
    }

    edge_tmp[which(phy_new$edge[,1] == idx1), 1] <- idx2
    edge_tmp[which(phy_new$edge[,2] == idx1), 2] <- idx2

    edge_tmp[which(phy_new$edge[,1] == idx2), 1] <- idx1
    edge_tmp[which(phy_new$edge[,2] == idx2), 2] <- idx1

    phy_new$edge <- edge_tmp
    return(phy_new)
}

demangle_ape <- function(phy) {
    phy_new <- reorder(phy, "postorder")
    ntip <- length(phy_new$tip.label)
    relabelled <- FALSE

    # Remove node labels
    if (is.null(phy_new$node.label))
    {
        warning("Removing Node Labels.")
        phy_new$node.label <- NULL
    }

    # Check for multifurcations produced by some phylogenetic reconstruction software and replace with binary mergers
    # In case this wasn't done already
    if (!is.binary(phy)) 
    {
        warning("Undated tree contains multifurcations. Mutlifurcations will be split into binary nodes prior to timing.")
        phy_new <- multi2di(phy_new)
    }

    # Generate new node labels
    phy_new <- makeNodeLabel(phy_new)
    phy_ape <- phy_new 
    if (is.rooted(phy_new))
    {
        # ape always assumes root has index ntip + 1. We certainly don't want this.
        root_pos_old <- ntip + 1
        # Move root to new position 
        root_pos_new <- 2*ntip - 1
        phy_new <- reindex_node(phy_new, root_pos_old, root_pos_new)
        relabelled <- TRUE
    }
    return(list(phy=phy_new, phy_ape=phy_ape, relabelled=relabelled))
}

unroot_demangled <- function(phyd) {
    # find branches incident to root
    phyd_u <- phyd
    br_arr <- cbind(phyd_u$edge, phyd_u$edge.length)
    ntip <- length(phyd_u$tip.label)
    r <- 2*ntip-1
    
    brs <- which(br_arr[,1] == r)

    br1 <- brs[1]
    br2 <- brs[2]

    len_br2 <- br_arr[br2,3]
    n2 <- br_arr[br2,2]

    br_arr[br1,1] <- n2 
    br_arr[br1,3] <- br_arr[br1,3] + len_br2

    br_root <- br_arr[br1,c(1,2)]
    br_arr <- br_arr[-br2, ]

    br_root_idx <- which((br_arr[,1] == br_root[1]) & (br_arr[,2] == br_root[2]))

    stopifnot(length(br_root_idx) == 1)

    phyd_u$edge <- br_arr[,c(1,2)]
    phyd_u$edge.length <- as.integer(br_arr[,3] + 0.5)

    len_merged <- br_arr[br_root_idx,3]

    return(list(unrooted=phyd_u, 
        root_br_pos=br_root_idx,
        split=(len_br2+0.5)/len_merged))
}

root_demangled <- function(phyd, root_br_pos, split) {
    # find branches incident to root
    phyd_r <- phyd
    br_arr <- cbind(phyd$edge, phyd$edge.length)
    ntip <- length(phyd$tip.label)
    r <- 2*ntip-1

    old_len <- br_arr[root_br_pos, 3]
    other_n <- br_arr[root_br_pos, 2]

    new_len <- floor(old_len * (1-split))
    resid_len <- old_len - new_len

    br_arr[root_br_pos, 2] <- r
    br_arr[root_br_pos, 3] <- new_len

    br_to_add <- cbind(other_n, r, resid_len)

    br_arr <- rbind(br_arr,br_to_add)
    phyd_r$edge <- br_arr[,c(1,2)]
    phyd_r$edge.length <- br_arr[,3]

    return(phyd_r)
}

mangle_ape <- function(phyd) {
    ntip <- length(phyd$tip.label)
    edges <- phyd$edge
    nset <- unique(as.vector(edges[,c(1,2)]))
    nset <- nset[nset>ntip]
    is_rooted <- min(sapply(nset, function(i) (length(which(edges[,1] == i )) + length(which(edges[,2] == i ))))) == 2

    if(!is.rooted)
    {
        warning("Output tree is not rooted!")
    }

    out <- phyd

    if(is_rooted) {
        r <- 2*ntip-1
        out <- reindex_node(phyd, ntip+1, r)
    }
    return(out)
}

dates_identical <- function(dates)
{
    return (var(dates, na.rm=T) == 0)
}

