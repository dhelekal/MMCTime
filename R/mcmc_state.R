init_tree_state <- function(phydata, branch_data, tip_dates_o, mu_init = 2, times_init=c())
{    
    times <- rep(0.0, phydata$n_node)
    qs <- as.integer(rep(0, phydata$n_node))
    taus <- rep(0.0, phydata$n_tip-1)
    n_tip <- phydata$n_tip

    times[1:n_tip] <- tip_dates_o
    tree_state <- list(taus=taus, times=times, qs=qs, n_tip=n_tip)
    if (length(times_init)!=0)
    {
        root_idx <- 2*n_tip-1
        ape_root <- n_tip+1
        tmp <- times_init[ape_root]
        times_init[ape_root] <- times_init[root_idx]
        times_init[root_idx] <- tmp 
        tree_state$times <- times_init
    } else 
    {
        initialise_times(phydata, tree_state, branch_data, mu_init)
    }
    taus_from_times(phydata, tree_state)
    return(tree_state)
}

init_par_state <- function(model)
{
    n <- model$n_par_tree + model$n_par_obs
    pars_constr <- rep(0.0, n)
    pars_unconstr <- runif(n, -2, 2)

    par_state <- list(pars=pars_unconstr, transf_pars=pars_constr)
    
    model$transform_pars(par_state)
    return(par_state)
}

get_par <- function(x, par_name, model)
{
    pnames <- c(model$names_par_obs, model$names_par_tree)
    idx <- which(pnames == par_name)
    if(length(idx) != 1)
    {
        stop(paste0("Parameter not found! Parameter names: ", paste0(pnames, collapse=" ")))
    }
    else
    {
        return(x$transf_pars[idx])
    }
}

make_trace <- function(model, draws, n_qs, n_taus, n_times, n_pa)
{
    names_taus <- paste0("tau_", 1:n_taus)
    names_times <- paste0("t_", 1:n_times)
    names_qs <- paste0("q_", 1:n_qs)
    names_pa <- paste0("pa_", 1:n_pa)
    
    names_pars <- c(model$names_par_obs, model$names_par_tree, model$names_summaries)
    cn <- c(names_taus, names_times, names_qs, names_pa, "root_br_pos", names_pars, "lp", ".iteration", ".chain")

    colnames(draws) <- cn
    return(draws)
}

post_process_trace <- function(trace, phydata)
{
    trace_cpy <- trace
    n_tip <- phydata$n_tip
    #Reindex root to ensure compatibility with ape!
    root_idx <- 2*n_tip - 1
    ape_root <- n_tip + 1
    n_node <- 2*n_tip-1

    vnames <- c("t_", "q_", "pa_")

    for (v in vnames)
    {
        trace_cpy[paste0(v, root_idx)] <- trace[paste0(v, ape_root)]
        trace_cpy[paste0(v, ape_root)] <- trace[paste0(v, root_idx)]
    }

    for (cn in paste0("pa_", 1:n_node))
    {
        trace_cpy[cn] <- sapply(trace_cpy[[cn]], function(i) if(!is.na(i) && i==ape_root) root_idx else if(!is.na(i) && i==root_idx) ape_root else i)
    }

    trace_cpy[paste0("tau_", root_idx-n_tip)] <- trace[paste0("tau_", ape_root-n_tip)]
    trace_cpy[paste0("tau_", ape_root-n_tip)] <- trace[paste0("tau_", root_idx-n_tip)]
    return(as_draws_df(trace_cpy))
}