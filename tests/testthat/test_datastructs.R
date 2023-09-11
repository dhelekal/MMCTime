context("Ape conversion")
test_that("Conversion between ape phylogenies and internal phylogeny works on a known tree",
{
    set.seed(2)

    n_samp <- c(2, 1, 1)
    samp_t <- c(0, 0.5, 0.53) 
    coal_t <- c(0.2, 0.52, 3.0)

    phy <- coal2ape(samp_t, n_samp, coal_t, c(2,2,2), node_name_prefix="N")$phy

    n <- length(n_samp)
    mu <- 100
    n_tip <- sum(n_samp)

    dates <- unlist(lapply(1:n, function(i) rep(samp_t[i], n_samp[i])))    

    phy$edge.length <- rep(5,length(phy$edge.length))
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata

    expect_equal(tr$sbounds, c(unname(dates[phy$tip.label]), 0, 0.5, 0.53))
    expect_equal(tr$n_tip, 4)
    expect_equal(tr$n_node, 7)

    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], mu)
    expect_equal(mcs$qs, rep(0,2*n_tip-1))
    expect_equal(mcs$times, c(unname(dates[phy$tip.label]),0.05, 0.55, 0.6))

    expect_true(validate_times(tr,mcs))
    expect_equal(tr, duplicate_tree(tr))
    
})

context("Parameter Handling")
test_that("Conversion between times and taus agrees",
{
    set.seed(2)

    n_samp <- c(2, 1, 1)
    samp_t <- c(0, 0.5, 0.53) 
    coal_t <- c(0.2, 0.52, 3.0)
    phy <- coal2ape(samp_t, n_samp, coal_t, c(2,2,2), node_name_prefix="N")$phy

    n <- length(n_samp)
    mu <- 100
    n_tip <- sum(n_samp)
    dates <- unlist(lapply(1:n, function(i) rep(samp_t[i], n_samp[i])))    
    phy$edge.length <- rep(5,length(phy$edge.length))
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], mu_init=mu)

    expect_true(validate_times(tr,mcs))

    taus_old <- c(mcs$taus)
    times_old <- c(mcs$times)

    taus_from_times(tr, mcs)
    taus_old <- c(mcs$taus)
    
    times_from_taus(tr, mcs)
    times_new <- c(mcs$times)

    taus_from_times(tr, mcs)
    taus_new <- c(mcs$taus)

    expect_equal(taus_new, taus_old)
    expect_equal(times_new, times_old)
})

context("Parameter Handling")
test_that("Deep copy functions work",
{
    set.seed(2)

    n_samp <- c(50, 50, 25, 25)
    samp_t <- c(0, 0.3, 0.4, 0.6) 
    
    sim_res <- sim_beta(samp_t, n_samp, 1.0, 0.2)
    phy <- coal2ape(samp_t, n_samp, sim_res$coal_times, sim_res$merger_sizes)$phy

    n <- length(n_samp)
    n_tip <- sum(n_samp)
    mu <- 50

    dates <- unlist(lapply(1:n, function(i) rep(samp_t[i], n_samp[i])))    
    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], mu_init=mu)
    taus_from_times(tr, mcs)

    mcs_copy <- duplicate_mc_state(mcs)
    expect_equal(mcs_copy, mcs)

    mcs_copy$qs[1] <- 1L
    expect_equal(mcs_copy$qs[1], 1)
    expect_equal(mcs$qs[1], 0)

    mc_state_copy_to(mcs, mcs_copy)

    expect_equal(mcs_copy, mcs)

    mcs_copy$qs[1] <- 1L
    expect_equal(mcs_copy$qs[1], 1)
    expect_equal(mcs$qs[1], 0)
})

context("Tree representation")
test_that("Internal binary tree representation agrees with binary tree simulator",
{
    set.seed(2)
    n_samp <- c(10, 15, 20)
    samp_t <- c(0, 0.3, 0.4) 

    mu <- 50

    samp_t_expanded <- expand_samp_times(samp_t, n_samp)
    samp_deltas <- rep(1, length(samp_t_expanded))

    n_samp_i <- as.integer(n_samp)
    sim_res <- sim_kingman(samp_t, n_samp_i, 1.0)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times

    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))

    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata
    mcs <- init_tree_state(tr, d$branch_lens, mu, times_init=phy_times)
    taus_from_times(tr, mcs)

    mcs_ts <-phy_times[order(phy_times)]# mcs$times[order(mcs$times)]

    evt_times_gt <- c(samp_t_expanded, sim_res$coal_times)
    delta_At_gt <- c(samp_deltas, -sim_res$merger_sizes+1)

    delta_At_gt <- delta_At_gt[order(evt_times_gt)]
    evt_times_gt <- evt_times_gt[order(evt_times_gt)]

    expect_equal(mcs_ts, evt_times_gt)

    as_coal_events <- as_lambda_events(tr, mcs)

    evt_times <- as_coal_events$event_times[order(as_coal_events$event_times)]
    delta_At <- as_coal_events$delta_At[order(as_coal_events$event_times)]

    expect_equal(evt_times, evt_times_gt)
    expect_equal(delta_At, delta_At_gt)
    expect_true(validate_times(tr,mcs))
})

test_that("Internal multi merger tree representation agrees with multi merger tree simulator",
{
    set.seed(3)
    n_samp <- c(100, 150, 50)
    samp_t <- c(0, 0.3, 0.4)

    mu <- 50

    samp_t_expanded <- expand_samp_times(samp_t, n_samp)
    samp_deltas <- rep(1, length(samp_t_expanded))

    n_samp_i <- as.integer(n_samp)
    
    sim_res <- sim_km_beta(samp_t, n_samp, 2/3, 2, 0.3)
    phy <- coal2ape(samp_t, n_samp_i, sim_res$coal_times, sim_res$merger_sizes)$phy
    
    phy <- multi2di(phy)
    phy_times <- node.depth.edgelength(phy)
    phy_times <- max(phy_times) - phy_times

    dates <- samp_t_expanded   
    names(dates) <- sapply(1:sum(n_samp), function(i) paste0("S", i))
    phy$edge.length <- rpois(length(phy$edge.length), mu * phy$edge.length)
    
    d <- tree_to_phydata(phy, dates)
    tr <- d$phydata
    mcs <- init_tree_state(tr, d$branch_lens, dates[phy$tip.label], mu_init=mu, times_init=phy_times)
    taus_from_times(tr, mcs)

    binary_to_mm(tr, mcs, tol=1e-8)
    taus_from_times(tr, mcs)

    evt_times_gt <- c(samp_t_expanded, sim_res$coal_times)
    delta_At_gt <- c(samp_deltas, -sim_res$merger_sizes+1)

    delta_At_gt <- delta_At_gt[order(evt_times_gt)]
    evt_times_gt <- evt_times_gt[order(evt_times_gt)]

    as_coal_events <- as_lambda_events(tr, mcs)

    evt_times <- as_coal_events$event_times[order(as_coal_events$event_times)]
    delta_At <- as_coal_events$delta_At[order(as_coal_events$event_times)]

    expect_equal(evt_times, evt_times_gt)
    expect_equal(delta_At, delta_At_gt)
    expect_true(validate_times(tr,mcs))
})