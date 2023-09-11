mcmc <- function (par_state_init, mcmc_state_init, phydata_init, branch_data, model, par_prior, mm_sd, enable_root_move)
{
    par_state_curr <- par_state_init
    mc_state_curr <- mcmc_state_init

    mc_state_prop <- duplicate_mc_state(mc_state_curr)
    par_state_prop <-  duplicate_par_state(par_state_curr)

    phydata_curr <- phydata_init
    phydata_prop <- duplicate_tree(phydata_curr)

    pivots <- find_pivots(phydata_curr, branch_data)

    n_tau <- length(mc_state_curr$taus)
    n_par <- length(par_state_curr$pars)

    tau_sd_scale <- .4
    par_sd_scale <- .20

    tau_sd_unsc <- rep(1/sqrt(n_tau), n_tau)
    par_sd_unsc <- rep(1/sqrt(n_par), n_par)

    tau_sd <- tau_sd_scale * tau_sd_unsc
    par_sd <- par_sd_scale * par_sd_unsc

    logJ_curr <- model$transform_pars(par_state_curr)
    logJ_prop <- model$transform_pars(par_state_prop)

    stopifnot(abs(logJ_curr-logJ_prop) < 1e-8)

    #lp <- function(mc_state_upd, par_state_upd, phy_upd) 
    #{
    #    lp_prior <- 0 
    #    lp_obs <- 0
    #    lp_prior <- model$prior_lp(if(phy_upd) phydata_prop else phydata_curr, if(mc_state_upd) mc_state_prop else mc_state_curr, if(par_state_upd) par_state_prop else par_state_curr)
    #    lp_obs <- model$obs_lp(if(phy_upd) phydata_prop else phydata_curr, if(mc_state_upd) mc_state_prop else mc_state_curr, if(par_state_upd) par_state_prop else par_state_curr, branch_data)
    #    return(lp_obs + lp_prior)
    #}

    lp2 <- function(mc_state_upd, par_state_upd, phy_upd) 
    {
        lp_pprior <- 0 
        lp_coal <- 0 
        lp_obs <- 0
        lp_pprior <- par_prior$prior_lp(if(par_state_upd) par_state_prop$transf_pars else par_state_curr$transf_pars) 
        lp_coal <- model$coal_lp(if(phy_upd) phydata_prop else phydata_curr, if(mc_state_upd) mc_state_prop else mc_state_curr, if(par_state_upd) par_state_prop else par_state_curr)
        lp_obs <- model$obs_lp(if(phy_upd) phydata_prop else phydata_curr, if(mc_state_upd) mc_state_prop else mc_state_curr, if(par_state_upd) par_state_prop else par_state_curr, branch_data)
        return(lp_pprior + lp_coal + lp_obs + if(par_state_upd) logJ_prop else logJ_curr)
    }
    
    #lp_curr <- lp(FALSE, FALSE, FALSE) 
    lp_curr <- lp2(FALSE, FALSE, FALSE) 

    lp_prop <- NA
    mh <- 0.0
    qr <- 0.0

    accept_mc <- function()
    {
        tmp <- mc_state_curr
        mc_state_curr <<- mc_state_prop
        mc_state_prop <<- tmp
        
        lp_curr <<- lp_prop
    }

    accept_par <- function()
    {
        tmp <- par_state_curr
        par_state_curr <<- par_state_prop
        par_state_prop <<- tmp
        
        logJ_curr <<- logJ_prop
        lp_curr <<- lp_prop
    }

    accept_tree <- function()
    {
        tmp <- phydata_curr
        phydata_curr <<- phydata_prop
        phydata_prop <<- tmp
        
        lp_curr <<- lp_prop
    }

    return(list
        (
            tau_move = function()
            {
                qr <<- rwm_tau_move(phydata_curr, mc_state_curr, mc_state_prop, tau_sd)
                lp_prop <<- lp2(TRUE, FALSE, FALSE)
                mh <<- lp_prop - lp_curr + qr
                if (log(runif(1)) < mh)
                {
                    accept_mc()
                    return(1)
                } else
                {
                    return(0)
                }
            },
            par_move = function()
            {
                qr <<- rwm_par_move(par_state_curr, par_state_prop, par_sd)
                logJ_prop <<- model$transform_pars(par_state_prop)

                lp_prop <<- lp2(FALSE, TRUE, FALSE)

                mh <<- lp_prop - lp_curr + qr
                if (log(runif(1)) < mh)
                {
                    accept_par()
                    return(1)
                } else
                {
                    return(0)
                }
            },
            mm_move = function()
            {

                qr <<- push_mm_move(phydata_curr, mc_state_curr, mc_state_prop, mm_sd, branch_data)
                lp_prop <<- lp2(TRUE, FALSE, FALSE)

                mh <<- lp_prop - lp_curr + qr
                if (log(runif(1)) < mh)
                {
                    accept_mc()
                    return(1)
                } else
                {
                    return(0)
                }
            },
            poly_move = function()
            {
                qr <<- topo_move(phydata_curr, mc_state_curr, phydata_prop, mc_state_prop, branch_data, pivots)

                lp_prop <<- lp2(TRUE, FALSE, TRUE)

                mh <<- lp_prop - lp_curr + qr
                if (log(runif(1)) < mh)
                {
                    accept_tree()
                    accept_mc()
                    return(1)
                } else
                {
                    return(0)
                }
            },
            root_move = ifelse(enable_root_move, function()
            {
                qr <<- root_move(phydata_curr, mc_state_curr, phydata_prop, mc_state_prop)

                lp_prop <<- lp2(TRUE, FALSE, TRUE)

                mh <<- lp_prop - lp_curr + qr
                if (log(runif(1)) < mh)
                {
                    accept_tree()
                    accept_mc()
                    pivots <<- find_pivots(phydata_curr, branch_data)
                    return(1)
                } else
                {
                    return(0)
                }
            }, function() stop("Root move not enabled!")),
            summaries = function()
            {
                model$summaries(phydata_curr, mc_state_curr, par_state_curr)
            },
            par_state = function()
            {
                par_state_curr
            },
            mc_state = function()
            {
                mc_state_curr
            },
            log_prob = function()
            {
                lp_curr
            },
            phydata = function()
            {
                phydata_curr
            },
            set_tau_scale = function(x)
            {
                tau_sd_scale <<- x
                tau_sd <<- tau_sd_unsc*tau_sd_scale
            },
            set_par_scale = function(x)
            {
                par_sd_scale <<- x
                par_sd <<- par_sd_unsc*par_sd_scale
            },
            get_tau_scale = function()
            {
                tau_sd_scale
            },
            get_par_scale = function()
            {
                par_sd_scale
            },
            set_tau_psd = function(x)
            {
                tau_sd_unsc <<- x
                tau_sd <<- tau_sd_unsc*tau_sd_scale
            },
            set_par_psd = function(x)
            {
                par_sd_unsc <<- x
                par_sd <<- par_sd_unsc*par_sd_scale
            },
            root_move_enabled = function()
            {
                enable_root_move
            }
        )
    )
}

warmup <- function(x, n_thin, n_eff_burnin, n_eff_var_est, n_eff_scale_est, chain_id=1, verbose=FALSE)
{
    burnin_it <- n_eff_burnin * n_thin
    var_est_it <- n_eff_var_est * n_thin
    sd_scale_it <- n_eff_scale_est * n_thin

    if(verbose) print(paste0("Chain ", chain_id, ": Starting warmup..."))
    poly_or_root <- NA
    if (x$root_move_enabled())
    {
        poly_or_root <- function() if(runif(1) <= 0.5) x$poly_move() else x$root_move()

    } else
    {
        poly_or_root <- function () x$poly_move()
    }
    for (i in 1:burnin_it)
    {
        poly_or_root()
        x$par_move()
        x$tau_move()
    }

    n_tau <- length(x$mc_state()$taus)
    n_par <- length(x$par_state()$pars)

    if(verbose) print(paste0("Chain ", chain_id, ": Estimating parameter variances..."))
    tau_mat <- matrix(rep(NA, n_tau * n_eff_var_est), nrow = n_tau, ncol = n_eff_var_est)
    k <- 1
    for (j in 1:var_est_it)
    {   
        poly_or_root()
        x$par_move()
        x$tau_move()

        if (j%%n_thin==0)
        {
            tau_mat[,k] <- x$mc_state()$taus
            k <- k+1
        }
    }
    x$set_tau_psd((apply(tau_mat, 1, sd)+0.05) / sqrt(n_tau))

    for (j in 1:burnin_it)
    {
        poly_or_root()
        x$par_move()
        if(runif(1) <= 0.5) x$tau_move() else x$mm_move()   
    }
    
    par_mat <- matrix(rep(NA, n_par * n_eff_var_est), nrow = n_par, ncol = n_eff_var_est)
    k <- 1
    for (j in 1:var_est_it)
    {
        poly_or_root()
        x$par_move()
        if(runif(1) <= 0.5) x$tau_move() else x$mm_move()


        if (j%%n_thin==0)
        {
            par_mat[,k] <- x$par_state()$pars
            k <- k+1
        }
    }
    x$set_par_psd((apply(par_mat, 1, sd)+0.05) / sqrt(n_par))
    
    if(verbose) print(paste0("Chain ", chain_id, ": Optimising proposal scalings..."))    
    
    n_acc_tau <- 0L
    n_acc_par <- 0L
    n_try_tau <- 0L

    target_accr <- 0.234
    ds <- 0.02

    for(j in 1:sd_scale_it)
    {
        poly_or_root()
        n_acc_par <- x$par_move() + n_acc_par
        if (runif(1) <= 0.5)
        {
            n_acc_tau <- x$tau_move() + n_acc_tau
            n_try_tau <- n_try_tau + 1L
        } else
        {
            x$mm_move()
        }
        if (j %% n_thin == 0)
        {
            accr_tau <- n_acc_tau / n_try_tau
            acrr_par <- n_acc_par / n_thin

            n_acc_tau <- 0L
            n_acc_par <- 0L

            n_try_tau <- 0L

            tsc <- x$get_tau_scale()

            if (accr_tau > target_accr)
            {
                tsc <- tsc + ds
            } else
            {
                tsc <- tsc - ds
            }

            stopifnot(tsc >= 0.05)
            x$set_tau_scale(tsc)

            psc <- x$get_par_scale()

            if (acrr_par > target_accr)
            {
                psc <- psc + ds
            } else
            {    
                psc <- psc - ds
            }
            
            stopifnot(psc >= 0.05)
            x$set_par_scale(psc)
        }
    }
} 

sample_one <- function(x, n_thin, chain_id, it, verbose=FALSE)
{
    moves <- c("par_move", "tau_move", "mm_move", "poly_move", "root_move")
    n_moves <- length(moves)
    n_acc <- rep(0L,n_moves)
    n_try <- rep(0L,n_moves)
    names(n_acc) <- moves
    names(n_try) <- moves

    do_move <- function(move)
    { 
        n_acc[[move]] <<- n_acc[[move]] + x[[move]]()
        n_try[[move]] <<- n_try[[move]] + 1L
    }
    if (x$root_move_enabled())
    {
        for (it in 1:n_thin)
        { 
            do_move(ifelse(runif(1) < 0.5, "poly_move", "root_move"))
            do_move("par_move")
            do_move(ifelse(runif(1) < 0.5, "tau_move", "mm_move"))
        }
    } else 
    {
        for (it in 1:n_thin)
        { 
            do_move("poly_move")
            do_move("par_move")
            do_move(ifelse(runif(1) < 0.5, "tau_move", "mm_move"))
        }
    }

    if (verbose)
    {
        accr <- n_acc/n_try
        print(paste0("Chain: ", chain_id," Acceptance rates: ", paste0(moves, ": ", accr, collapse=" ")))
    }
    
    mcs <- x$mc_state()
    tre <- x$phydata()
    ps <- x$par_state()
    sums <- x$summaries()

    c(list(), mcs$taus, mcs$times, mcs$qs, tre$topo_mat[,1], tre$root_br_pos, ps$transf_pars, sums, x$log_prob(), it, chain_id)
}