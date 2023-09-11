tree_to_phydata <- function(phy, dates, enable_root_move=T)
{
    # Relabel tree so that root can be at positon 2n-1
    n_tip <- length(phy$tip.label)
    n_node <- 2*n_tip - 1
    
    # Required for obs negbin distributions
    phy$edge.length <- round(phy$edge.length)

    dem <- demangle_ape(phy)
    phy_ape <- dem$phy_ape
    
    if (is.rooted(dem$phy) && enable_root_move) 
    {
        ur <- unroot_demangled(dem$phy)
        phy_d <- ur$unrooted
        root_br_pos <- as.integer(ur$root_br_pos)
    } else if(is.rooted(dem$phy))
    {
        phy_d <- dem$phy
        root_br_pos <- as.integer(NA)
    } else
    {
        stop("Tree must be rooted. Consider using outgroup rooting to root the tree.")
    }

    br_mat <- cbind(phy_d$edge, phy_d$edge.length)
    colnames(br_mat) <- c("e1", "e2", "mut")
    mode(br_mat) <- "integer"
    
    dates_o <- dates[phy$tip.label]
    internal_set <- (n_tip+1):n_node
    tip_set <- 1:n_tip

    topo_mat <- build_topo_mat(br_mat, n_node, root_br_pos, 1e-6)
    sbounds <- find_sbounds(topo_mat, dates_o, n_tip, n_node-1)

    #return(list(edges=br_mat, 
    #    incident_nodes=incident_nodes,
    #    sbounds=sbounds,
    #    parent_branches=parent_branches,
    #    root_br_pos=root_br_pos,
    #    n_tip=n_tip,
    #s    n_node=n_node))

    phydata <- list(topo_mat=topo_mat, 
        sbounds=sbounds,
        root_br_pos=root_br_pos,
        n_tip=n_tip,
        n_node=n_node,
        root_fixed=ifelse(enable_root_move, 0L, 1L))
    branch_lens <- br_mat[,3]

    return(list(phydata=phydata,
        branch_lens=branch_lens,
        phy_ape=unroot(phy_ape)
    ))
}   

validate_dates <- function(phy, dates)
{
    stopifnot("Sampling dates cannot have missing names!"=all(!is.null(names(dates)) & !is.na(names(dates))))
    stopifnot("Sampling dates names must all be unique!"=length(names(dates)) == length(unique(names(dates))))
    stopifnot("Length of sampling date names must match number of tips!"=length(names(dates)) == length(phy$tip.label))

    for (dn in names(dates))
    {
        if (!(dn %in% phy$tip.label)) stop(paste0("Cannot find a matching tip label for date name: ", dn))
    }
}