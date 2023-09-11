context("Coalescent likelihoods")
test_that("Kingman likelihood matches that returned by simulator",
{
    set.seed(4)
    n_samp <- c(10, 15, 20)
    samp_t <- c(0, 0.3, 0.4) 

    mu <- 50
    samp_t_expanded <- expand_samp_times(samp_t, n_samp)

    nu <- 0.14

    n_samp_i <- as.integer(n_samp)
    sim_res <- sim_kingman(samp_t, n_samp_i, nu)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times

    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata    
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], times_init=phy_times)

    lp <- kingman_lp(tr,mcs,nu)

    expect_equal(sim_res$sim_lp, lp)    
})

test_that("KM-BETA likelihood matches that returned by simulator",
{
    set.seed(4)
    n_samp <- c(50, 80, 90)
    samp_t <- c(0, 0.3, 0.8) 

    mu <- 50
    phi <- 3/4
    nu <- 4
    alpha <- 0.3
    samp_t_expanded <- expand_samp_times(samp_t, n_samp)

    n_samp_i <- as.integer(n_samp)
    sim_res <- sim_km_beta(samp_t, n_samp, phi, nu, alpha)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    phy <- multi2di(phy)
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times

    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata    
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], times_init=phy_times)
    
    binary_to_mm(tr, mcs, tol=1e-8)
    taus_from_times(tr, mcs)
    
    kmb_lp <- km_beta_lp(tr, mcs, nu, phi, alpha)

    expect_equal(sim_res$sim_lp, kmb_lp)   
})

test_that("BETA likelihood matches that returned by simulator",
{
    set.seed(4)
    n_samp <- c(50, 80, 90)
    samp_t <- c(0, 0.3, 0.8) 

    mu <- 50
    nu <- 4
    alpha <- 0.3
    samp_t_expanded <- expand_samp_times(samp_t, n_samp)

    n_samp_i <- as.integer(n_samp)
    sim_res <- sim_beta(samp_t, n_samp, nu, alpha)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    phy <- multi2di(phy)
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times

    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata    
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], times_init=phy_times)
    
    binary_to_mm(tr, mcs, tol=1e-8)
    taus_from_times(tr, mcs)
    
    beta_lp <- beta_lp(tr, mcs, nu, alpha)

    expect_equal(sim_res$sim_lp, beta_lp)   
})

test_that("DS likelihood matches that returned by simulator",
{
    set.seed(4)
    n_samp <- c(50, 80, 90)
    samp_t <- c(0, 0.3, 0.8) 

    mu <- 50
    nu <- 4
    phi <- .8
    samp_t_expanded <- expand_samp_times(samp_t, n_samp)

    n_samp_i <- as.integer(n_samp)
    sim_res <- sim_ds(samp_t, n_samp, phi, nu)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    phy <- multi2di(phy)
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times
    
    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata    
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], times_init=phy_times)
    
    binary_to_mm(tr, mcs, tol=1e-12)
    taus_from_times(tr, mcs)
    
    ds_lp <- ds_lp(tr, mcs, nu, phi)

    expect_equal(sim_res$sim_lp, ds_lp)   
})

context("Parameter transformations")
test_that("Jacobian for tau transform and it's inverse cancel",
{
    set.seed(3)
    n_samp <- c(10, 15, 20)
    samp_t <- c(0, 0.3, 0.4) 

    mu <- 50

    samp_t_expanded <- expand_samp_times(samp_t, n_samp)
    samp_deltas <- rep(1, length(samp_t_expanded))

    n_samp_i <- as.integer(n_samp)
    
    sim_res <- sim_beta(samp_t, n_samp_i, 1.0, .2)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    
    phy <- multi2di(phy)
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times
    
    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))
    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata    
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], times_init=phy_times)
        
    log_j <- taus_to_times_logJ(tr, mcs)
    log_j_inv <- times_to_taus_logJ(tr, mcs)

    expect_equal(log_j, -log_j_inv)

    binary_to_mm(tr, mcs, tol=1e-8)
    log_j_mm <- taus_to_times_logJ(tr, mcs)
    log_j_inv_mm <- times_to_taus_logJ(tr, mcs)

    expect_equal(log_j_mm, -log_j_inv_mm)
})


context("Branch likelihoods")
test_that("Computed branch likelihood matches simulator likelihood",
{
    set.seed(2)

    mu <- 10.2
    omega<-3.3

    n_samp <- c(2, 1, 1)
    samp_t <- c(0, 0.5, 1.0) 
    coal_t <- c(0.25, 0.75, 2.0)

    phy <- coal2ape(samp_t, n_samp, coal_t, c(2,2,2), node_name_prefix="N")$phy
    phy_mut <- phy
    phy_mut$edge.length <- rep(5,length(phy_mut$edge.length))

    br_root_1 <- 1.0
    br_root_2 <- 1.25
    br_oth <- c(0.25,0.25,0.50,0.25)

    nb_lp <- function(x,v,mu,omega,l) 
    {
        theta <- mu / omega
        p <- 1/(1+omega)
        dnbinom(x, theta * v, rep(p,length(v)), log=l)
    }


    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times

    dates <- expand_samp_times(samp_t, n_samp)   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    for (m in MODEL_OPTS())
    {
        mod<-tree_model(m)
        ps <- init_par_state(mod)

        mu <- get_par(ps, "mu", mod)
        omega <- get_par(ps, "omega", mod)

        lp_gt_ur <- sum(nb_lp(rep(5,4), br_oth, mu, omega, T)) + log(sum(sapply(0:10, function(i) nb_lp(10-i, br_root_2, mu, omega, F) * nb_lp(i, br_root_1, mu, omega, F))))
        lp_gt_r <- sum(nb_lp(rep(5,6), c(br_root_1, br_root_2, br_oth), mu, omega, T))

        d_ur <- tree_to_phydata(phy_mut, dates)
        mcs_ur <- init_tree_state(d_ur$phydata, d_ur$branch_lens, dates[phy$tip.label], times_init=phy_times)
        lp_comp_ur <- mod$obs_lp(d_ur$phydata, mcs_ur, ps, d_ur$branch_lens)
        expect_equal(lp_comp_ur,lp_gt_ur)

        d_r <- tree_to_phydata(phy_mut, dates, F)
        mcs_r <- init_tree_state(d_r$phydata, d_r$branch_lens, dates[phy$tip.label], times_init=phy_times)
        lp_comp_r <- mod$obs_lp(d_r$phydata, mcs_r, ps, d_r$branch_lens)
        expect_equal(lp_comp_r,lp_gt_r)
    }
})