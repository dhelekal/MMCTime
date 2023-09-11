#' Time an ML tree under a Lambda (multiple-merger) coalescent
#' @param phy A rooted ML tree, with branch lengths scaled in the absolute number of substitutions (not per site!). Rooting can be arbitrary.
#' @param dates A vector of tip dates, must share names with tip labels
#' @param n_draws (Optional) Number of draws per chain to retain. Default: 1e3
#' @param n_thin (Optional) Number of iterations per retained draw. Default: 1e3
#' @param n_chain (Optional) Number of chains to run in parallel. Default: 4
#' @param mm_sd (Optional) Multiple merger move sigma. Default: 0.7
#' @param model Lambda Coalescent prior selection. One of `beta` - Beta Coalescent; `km_beta` - Extended Beta Coalescent; `kingman` - Kingman's Coalescent (Binary)
#' @param verbose (Optional) verbose. Default: FALSE
#' @param fix_root (Optional) is the root position assumed to be known. Default: FALSE
#' @param ... (Optional) Additional arguments passed to parameter priors. See `PPAR_NAMES` for available prior hyperparameters and `PPAR_DEFAULTS` for default values.
#' @export
mmctime <- function(phy, dates, n_draws=1e3, thin=1e3, n_chain=4, mm_sd=0.7, model, verbose=FALSE, fix_root=FALSE, ...) 
{
    # Sanity checks
    stopifnot("Dates must be numeric"=is.numeric(dates))
    stopifnot("length(dates) must match the number of tips"=length(dates)==length(phy$tip.label))
    stopifnot("Branch lengths must be in number of mutations per genome, not per site!" =max(phy$edge.length) >= 1)

    N <- length(phy$tip.label)
    if (thin < N)
    {
        warning("Number of iterations per draw is lower than the number of internal nodes in phylogeny. This can lead to high autocorrelations.")
    }

    n_iter <- as.integer(n_draws * thin)
    stopifnot("Number of iterations must be greater than 0" = n_iter > 0)

    # Convert to coalescent time
    dmax <- max(dates)
    dates_ct <- dmax - dates
    
    # Make sure that the labels are consistent with the phylogeny
    validate_dates(phy, dates_ct)

    # Reorder dates to match tip label order
    dates_ct <- dates_ct[phy$tip.label]

    # Convert ape phylogeny to internal representation
    phy_conv <- tree_to_phydata(phy, dates_ct, !fix_root)
    phydata <- phy_conv$phydata
    branch_lens <- phy_conv$branch_lens

    undated <- phy_conv$phy_ape

    #Generate model object
    mod <- tree_model(model)

    #Generate parameter priors
    pnames <- c(mod$names_par_obs, mod$names_par_tree)
    par_prior <- make_par_prior(pnames, model, ...)

    #Chain initialisation function 
    init_chain <- function()
    {
        par_state <- init_par_state(mod)   
        mc_state <- init_tree_state(phydata, branch_lens, dates_ct, mu_init = par_state$transf_pars[1]) 
        chain <- mcmc(par_state, mc_state, phydata, branch_lens, mod, par_prior, mm_sd, !fix_root)
        chain
    }
    
    draws <- NA
    par_seed <- sample.int(10^8,1)
    if (n_chain>1)
    {
        cores <- min(detectCores(logical = TRUE) - 1, n_chain) 
        cl <- makeCluster(cores, outfile="")
        registerDoParallel(cl)
        registerDoRNG(seed = par_seed)

        draws <- foreach(i=1:n_chain, .combine = 'c') %dopar% 
        {
            chain <- init_chain()
            warmup(chain, thin, 100, 100, 50, chain_id=i, verbose=verbose)
            ld <- lapply(1:n_draws, function(it) sample_one(chain, thin, i, thin*(it-1L)+1L, verbose=verbose))
            
            tmp_df <- do.call(rbind.data.frame, ld)
            colnames(tmp_df)<-1:ncol(tmp_df)
            list(tmp_df)
        }      
        stopCluster(cl)
        draws<-do.call(rbind, draws)
    } else
    {
        chain <- init_chain()
        warmup(chain, thin, 100, 100, 50, chain_id=1L, verbose=verbose)
        ld <- lapply(1:n_draws, function(it) sample_one(chain, thin, 1L, thin*(it-1L)+1L, verbose=verbose))
        draws <- do.call(rbind.data.frame, ld)
    }

    n_qs <- 2*N-1
    n_times <- 2*N-1
    n_pa <- 2*N-1
    n_taus <- N-1

    trace <- make_trace(mod, draws, n_qs, n_taus, n_times, n_pa)
    return(timingRes(post_process_trace(trace, phydata), undated, mod$names_par_obs, mod$names_par_tree, mod$names_summaries, par_prior))
}
