timingRes <- function(draws_df, undated, sampling_dates, names_par_obs, names_par_tree, names_summaries, par_prior)
{
    n_draws <- nrow(draws_df)
    root_n <- length(undated$tip.label) + 1
    draws_df$tree_length <- sapply(1:n_draws, function(i) comp_treelen(sample2tree_internal(draws_df, undated, i)))
    s <- summarise_draws(draws_df)
    s <- s[s$variable %in% c(paste0("t_", root_n), names_par_obs, names_par_tree, names_summaries, "tree_length"), ]

    out <- list(draws=draws_df, 
        undated=undated, 
        sampling_dates=sampling_dates,
        n_draws=n_draws,
        summaries=s, 
        names_par_obs=names_par_obs, 
        names_par_tree=names_par_tree, 
        names_summaries=names_summaries,
        par_prior=par_prior)
    attr(out, "class") <- "timingRes"
    return(out)
}

#' @export
sample_timetree <- function(x, ...) UseMethod("sample_timetree")
#' @export 
plot_treeCI <- function(x, ...) UseMethod("plot_treeCI")
#' @export 
plot_traces <- function(x, ...) UseMethod("plot_traces")
#' @export 
plot_trace <- function(x, ...) UseMethod("plot_trace")
#' @export 
plot_pars <- function(x, ...) UseMethod("plot_pars")
#' @export 
plot_densiCI <- function(x, ...) UseMethod("plot_densiCI")
#' @export
plot_mm_tree <- function(x, ...) UseMethod("plot_mm_tree")
#' @export
build_cons_tree <- function(x, ...) UseMethod("build_cons_tree")

#' @export
summary.timingRes <- function(x)
{
    cat("\nMutlimerger tree timing result:", deparse(substitute(x)), "\n\n")
    cat("Summary:\n")
    print(x$summaries)
}

#' @export
print.timingRes <- function(x)
{
    cat("\nMutlimerger tree timing result:", deparse(substitute(x)), "\n\n")
    cat("Summary:\n")
    print(x$summaries)
}

#' Sample time tree from the posterior
#' @param x an object of class timingRes
#' @param n_samp number of samples to draw
#' @param replace sample with replacement
#' @return an object of class `ape::phylo` or `ape::multiPhylo`
#' @export
sample_timetree.timingRes  <- function(x, n_samp=1, replace=TRUE)
{
    stopifnot(replace || (n_samp <= x$n_draws))
    draws <- sample.int(x$n_draws, n_samp,replace=replace)
    out <- NA
    if (n_samp == 1)
    {
        out <- di2multi(sample2tree(x, draws))
        out$root.edge <- 0
    } else
    {
        out <- lapply(draws, function(j) {
            t <- di2multi(sample2tree(x, j))
            t$root.edge <- 0
            t
        })
        class(out) <- "multiPhylo"
    }
    
    return(out)
}

#' Build a majority rule consensus tree
#' @param x an object of class timingRes
#' @param compute_summaries should clade height quantiles be computed
#' @return a list of: `phy` - the majority rule tree, `support` - clade supports, `ts_qs` - clade height quantiles 
#' @export
build_cons_tree.timingRes <- function(x, compute_summaries = FALSE)
{
    tres_bi <- lapply(1:x$n_draws, function(i) sample2tree(x, i))
    tres_mm <- lapply(tres_bi, di2multi)

    class(tres_bi) <- "multiPhylo"
    class(tres_mm) <- "multiPhylo"

    edge <- maj_edge(tres_mm)
    tiplabs <- x$undated$tip.label
    n_edge <- nrow(edge)
    n_tip <- length(tiplabs)
    n <- n_edge + 1
    m <- 2*n_tip - 1
    elen_tmp <- rep(1, n_edge) 

    tr <- list(edge=edge, edge.length=elen_tmp, Nnode=n_edge - n_tip + 1, tip.label=tiplabs)
    class(tr) <- "phylo" 
    
    clabs <- lapply(1:n, function(i) if(i > n_tip) extract.clade(tr, i)$tip.label else tr$tip.label[i])
    clade_idx <- find_clade_index(tres_bi, clabs)

    t_draws <- as.matrix(suppressWarnings(x$draws[, paste0("t_", 1:m)]))
    clade_ts <- matrix(NA, x$n_draws, n)
    for (i in 1:x$n_draws) 
    {
        clade_ts[i,] <- t_draws[i,clade_idx[i,]]
    }
    
    med_ts <- apply(clade_ts,2,median, na.rm=T)
    support <- apply(!is.na(clade_idx),2,sum)/x$n_draws

    tr <- times_to_elen(tr, med_ts)
    if(!compute_summaries)
    {
        return(list(phy=tr, support=support))
    } else
    {
        ts_qs <- t(apply(clade_ts, 2, function(y) quantile(y, c(0.025, .5, 0.975), na.rm=T)))
        ts_qs[1:n_tip, c(1,3)] <- NA
        return(list(phy=tr, support=support, ts_qs=ts_qs))
    }
}

#' Plot a modified densiTree
#' @param x an object of class timingRes
#' @param mrsd (Optional) most recent sampling time. Default: determine from tip times.
#' @param n_samp (Optional) number of draws to overlay in the densiTree. Default: 100
#' @param layout (Optional) tree layout passed to ggtree. Default: "rectangular"
#' @param tip.order (Optional) a vector specifying the order in which the tips should be displayed. Default: determine from data
#' @export
plot_densiCI.timingRes <- function(x, mrsd=NULL, n_samp=100, layout="rectangular", tip.order=NULL)
{
    densiCI(x, mrsd, n_samp, layout, tip.order)
}

#' Plot parameter pair marginals
#' @param x an object of class timingRes
#' @export
plot_pars.timingRes <- function(x)
{
    root_n <- length(x$undated$tip.label) + 1
    pnames<-c(x$names_par_obs, paste0("t_", root_n), x$names_par_tree)
    mcmc_pairs(subset_draws(x$draws,variable=pnames), pars=pnames, off_diag_fun = "hex") 
}

#' Plot parameter traces
#' @param x an object of class timingRes
#' @export
plot_traces.timingRes <- function(x)
{
    root_n <- length(x$undated$tip.label) + 1
    pnames<-c(x$names_par_obs, paste0("t_", root_n), x$names_par_tree, x$names_summaries)
    mcmc_trace(subset_draws(x$draws,variable=pnames), pars=pnames, facet_args = list(ncol = 1, strip.position = "left")) 
}

#' Plot a multiple merger tree
#' @param x an object of class `ape::phylo`
#' @export
plot_mm_tree.phylo <- function(x, ...)
{
    x <- di2multi(x)
    
    n <- length(x$tip.label)
    m <- nrow(x$edge) + 1

    msize_node <- unname(c(rep(1,n),table(x$edge[,1])))
    msize_outgoing <- rep(F, m)

    for (i in 1:(m-1))
    {
        msize_outgoing[x$edge[i,2]] <- msize_node[x$edge[i,1]]
    }
    
    phy <- x %>% as_tibble %>% mutate(
        msize = msize_outgoing
    ) %>% as.treedata

    plt <- ggtree(phy, aes(color=msize>2), ...)
    plt <- plt +
        scale_color_manual(values=c("red","blue")) +
        labs(color="Multiple Merger") +
        theme_tree2()
    plt
}

#' Plot several multiple merger trees
#' @param x an object of class `ape::multiPhylo`
#' @export
plot_mm_tree.multiPhylo <- function(x, ...)
{
    plt <- wrap_plots(lapply(x, function(xx) plot_mm_tree(xx, ...))) + plot_layout(guides = 'collect')
    plt
}