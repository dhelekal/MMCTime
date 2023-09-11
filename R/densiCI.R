## ggtree internal function
## All credit to Guangchuang Yu
getYcoord <- function(tr, step=1, tip.order = NULL) {
    Ntip <- length(tr[["tip.label"]])
    N <- getNodeNum(tr)

    edge <- tr[["edge"]]
    parent <- edge[,1]
    child <- edge[,2]

    cl <- split(child, parent)
    child_list <- list()
    child_list[as.numeric(names(cl))] <- cl

    y <- numeric(N)
    if (is.null(tip.order)) {
        tip.idx <- child[child <= Ntip]
        y[tip.idx] <- 1:Ntip * step
    } else {
        tip.idx <- 1:Ntip
        y[tip.idx] <- match(tr$tip.label, tip.order) * step
    }
    y[-tip.idx] <- NA


    pvec <- edge2vec(tr)

    currentNode <- 1:Ntip
    while(anyNA(y)) {
        ## pNode <- unique(parent[child %in% currentNode])
        pNode <- unique(pvec[currentNode])

        ## piping of magrittr is slower than nested function call.
        ## pipeR is fastest, may consider to use pipeR
        ##
        ## child %in% currentNode %>% which %>% parent[.] %>% unique
        ## idx <- sapply(pNode, function(i) all(child[parent == i] %in% currentNode))
        idx <- sapply(pNode, function(i) all(child_list[[i]] %in% currentNode))
        newNode <- pNode[idx]

        y[newNode] <- sapply(newNode, function(i) {
            mean(y[child_list[[i]]], na.rm=TRUE)
            ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
        })

        currentNode <- c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
        ## currentNode <- c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
        ## parent %in% newNode %>% child[.] %>%
        ##     `%in%`(currentNode, .) %>% `!` %>%
        ##         currentNode[.] %>% c(., newNode)
    }

    return(y)
}

## ggtree internal function
## All credit to Guangchuang Yu
edge2vec <- function(tr) {
  parent <- tr$edge[,1]
  child <- tr$edge[,2]
  
  ## use lookup table
  pvec <- integer(max(tr$edge))
  pvec[child] <- parent
  return(pvec)
}

## densi_CI internal function
## densitree but for majority rule consensus clades we:
## use median heights, plot credible intervals and align, y positions
densiCI <- function(x, mrsd, n_samp, layout, tip.order)
{
    trees_bi <- lapply(1:x$n_draws, function(i) sample2tree(x, i))

    cons <- build_cons_tree(x, T)
    maj_tree <- cons$phy

    tiplabs <- x$undated$tip.label
    n_tip <- length(tiplabs)

    n <- nrow(maj_tree$edge)
    m <- (2*n_tip-1)

    if(is.null(mrsd)) mrsd <- max(x$sampling_dates)

    qs <- cons$ts_qs
    colnames(qs) <- c("lower", "median", "upper")

    support <- cons$support
    support <- round(support, 2)
    support[support < .5] <- NA

    clabs <- lapply(1:n, function(i) if(i > n_tip) extract.clade(maj_tree, i)$tip.label else maj_tree$tip.label[i])
    c_med_times <- cons$ts_qs[,2]
    offset <- min(c_med_times)

    trees <- sample_timetree(x, n_samp)
    
    m_clade_idx <- find_clade_index(trees, clabs)


    # Adjust heights of majority clades to match median
    # Visualise uncertainty using CIs for these later
    # This makes the plot legible
    for (i in 1:n_samp)
    {
        nh <- node.depth.edgelength(trees[[i]])
        nh <- max(nh)-nh + offset
        for(j in 1:n)
        {
            if (!is.na(m_clade_idx[i,j]))
            {
                nh[m_clade_idx[i,j]] <- c_med_times[j]
            }
        }
        trees[[i]]$edge.length <- sapply(1:nrow(trees[[i]]$edge), function(j) nh[trees[[i]]$edge[j,1]] - nh[trees[[i]]$edge[j,2]])
    }


    ### Begin code from ggtree/ggdensitree.R
    ### All credit to Guangchuang Yu
	trees <- lapply(trees, as.phylo)
	trees.f <- suppressWarnings(lapply(trees, fortify, layout=layout))

    if (is.null(tip.order))
    {
        dist_f <- function(tr)
        {
            ndepths <- node.depth(tr, 1)
            m <- mrca(tr)
            for (i in 1:nrow(m))
            {
                m[i,] <- ndepths[m[i,]]
            }
            as.dist(m-1)
        }

        dists <- lapply(trees, dist_f)
        h <- hclust(Reduce("+", dists)/n_samp)
        tip.order <- h$label[h$order]

    }
    ### End code from ggdensitree.R

    maj_tree.f <- fortify(as.phylo(maj_tree), layout=layout)

    max_x_med <- max(maj_tree.f$x,na.rm=TRUE)
    maj_tree.f$y <- getYcoord(maj_tree, tip.order = tip.order)
    maj_tree.f$x <- maj_tree.f$x - max_x_med

    # Add CIs for majority rule tree
    med_tr <- maj_tree.f %>% as_tibble %>% mutate(time_med = as.list(cons$ts_qs[,2]), 
        time_range = apply(qs[,c(1,3)],1,as.list),
        support = as.list(support)
    )
    
    ### Begin modified code from ggtree/ggdensitree.R
    max.x <- vapply(trees.f, function(y) max(y$x, na.rm = TRUE), numeric(1))
	trees.f <- lapply(1:length(trees), function(i) {
		trees.f[[i]]$y <- getYcoord(trees[[i]], tip.order = tip.order)
		trees.f[[i]]$x <- trees.f[[i]]$x - max.x[i]
        
        ## Align majority clade vertical positions
        for (j in 1:length(m_clade_idx[i,]))
        {
            c_idx <- m_clade_idx[i,j]
            if(!is.na(c_idx)) trees.f[[i]]$y[c_idx] <- maj_tree.f$y[j]
        }
        
		trees.f[[i]]
	})  
    ### End modified code from ggdensitree.R
    
    trees.f <- lapply(1:length(trees), function(i) 
    {
        n_edge <- nrow(trees[[i]]$edge)
        msize_node <- unname(c(rep(1,n_tip),table(trees[[i]]$edge[,1])))
        msize_outgoing <- rep(F, n_edge + 1)
        msize_outgoing[trees[[i]]$edge[1:n_edge,2]] <- msize_node[trees[[i]]$edge[1:n_edge,1]]
        tr <- trees.f[[i]] %>% as_tibble %>% mutate(msize = msize_outgoing)
        tr
    })

    # Plot CIs for majority rule tree
    # Keep the tree hidden except for node positions and CIs
    plt <- ggtree(tr=med_tr, layout=layout, color=NA, alpha=0.0)
    # Overlay densitrees

    plt <- plt + 
        lapply(trees.f, function(tr) geom_tree(data=tr, layout=layout, aes(color=msize>2), alpha=1/n_samp)) |> ggblend::blend("saturate") 

    plt <- plt + geom_range(range="time_range", center="time_med", color="gray65", size = 1, alpha = 0.8) +
            geom_nodepoint(size = .4, color="black") +
            geom_tippoint(size = .25, color="black") +
            geom_nodelab(aes(label=support), vjust=-.5, size=1.5)

    plt + 
    theme_tree2() + 
    scale_color_manual(values=c("red","blue")) + 
    labs(color="Multiple Merger") +
    scale_x_continuous(labels = scales::label_math(expr=.x, format=function(x) mrsd + x)) +
    theme(
        axis.text.x=element_text(angle = 45, hjust=1, size=4),
        panel.grid.major = element_blank(), 
        axis.line = element_line(size=rel(0.2), colour = "grey80"),
        plot.title = element_text(hjust = 0.5,size=6.0)
    )
}