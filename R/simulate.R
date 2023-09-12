#' Simulate a phylogeny from the extended beta coalescent
#' @param samp_times sampling epochs
#' @param n_samp number of samples per epoch
#' @param phi phi
#' @param nu nu
#' @param alpha alpha^*
#' @return a list consisting of sampling dates for each tip and the resulting phylogeny
#' @export
simulate_km_beta <- function(samp_times, n_samp, phi, nu, alpha)
{
    stopifnot("alpha must be between 0 and 1"= alpha<=1 && alpha > 0)
    stopifnot("phi must be between 0 and 1"= phi<=1 && phi >= 0)
    stopifnot("nu must be positive"= nu > 0)

    n_samp_i <- as.integer(n_samp)

    tmax <- max(samp_times)
    samp_times <- tmax - samp_times    
    n_samp <- n_samp[order(samp_times)]
    samp_times <- samp_times[order(samp_times)]

    validate_input(samp_times, n_samp_i)
    res <- sim_km_beta(samp_times, n_samp_i, phi, nu, alpha)

    tmp <- coal2ape(samp_times, n_samp_i, res$coal_times, res$merger_sizes)
    phy <- tmp$phy 
    dates <- tmax - tmp$dates

    return(list(phy=phy, dates=dates, lp=res$sim_lp))
}
#' Simulate a phylogeny from a particular durrett and schweinsberg coalescent. Experimental.
#' @param samp_times sampling epochs
#' @param n_samp number of samples per epoch
#' @param phi phi
#' @param nu nu
#' @return a list consisting of sampling dates for each tip and the resulting phylogeny
#' @export
simulate_durret_schweinsberg <- function(samp_times, n_samp, phi, nu)
{
    stopifnot("phi must be between 0 and 1"= phi<=1 && phi >= 0)
    stopifnot("nu must be positive"= nu > 0)

    n_samp_i <- as.integer(n_samp)

    tmax <- max(samp_times)
    samp_times <- tmax - samp_times    
    n_samp <- n_samp[order(samp_times)]
    samp_times <- samp_times[order(samp_times)]

    validate_input(samp_times, n_samp_i)
    res <- sim_ds(samp_times, n_samp_i, phi, nu)

    tmp <- coal2ape(samp_times, n_samp_i, res$coal_times, res$merger_sizes)
    phy <- tmp$phy 
    dates <- tmax - tmp$dates

    return(list(phy=phy, dates=dates, lp=res$sim_lp))
}

#' Simulate a phylogeny from the beta coalescent
#' @param nu nu
#' @param alpha alpha^*
#' @return a list consisting of sampling dates for each tip and the resulting phylogeny
#' @export
simulate_beta <- function(samp_times, n_samp, nu, alpha)
{
    stopifnot("nu must be positive"= nu > 0)
    stopifnot("alpha must be between 0 and 1"= alpha<=1 && alpha > 0)
    n_samp_i <- as.integer(n_samp)

    tmax <- max(samp_times)
    samp_times <- tmax - samp_times    
    n_samp <- n_samp[order(samp_times)]
    samp_times <- samp_times[order(samp_times)]

    validate_input(samp_times, n_samp_i)
    res <- sim_beta(samp_times, n_samp_i, nu, alpha)

    tmp <- coal2ape(samp_times, n_samp_i, res$coal_times, res$merger_sizes)
    phy <- tmp$phy 
    dates <- tmax - tmp$dates

    return(list(phy=phy, dates=dates, lp=res$sim_lp))
}

#' Simulate a phylogeny from kingmans coalescent
#' @param samp_times sampling epochs
#' @param n_samp number of samples per epoch
#' @param nu nu
#' @return a list consisting of sampling dates for each tip and the resulting phylogeny
simulate_kingman <- function(samp_times, n_samp, nu)
{
    stopifnot("nu must be positive"= nu > 0)
    n_samp_i <- as.integer(n_samp)

    tmax <- max(samp_times)
    samp_times <- tmax - samp_times    
    n_samp <- n_samp[order(samp_times)]
    samp_times <- samp_times[order(samp_times)]

    validate_input(samp_times, n_samp_i)
    res <- sim_kingman(samp_times, n_samp_i, nu)

    tmp <- coal2ape(samp_times, n_samp_i, res$coal_times, res$merger_sizes)
    phy <- tmp$phy 
    dates <- tmax - tmp$dates

    return(list(phy=phy, dates=dates, lp=res$sim_lp))
}

#' Simulate a mutation scaled phylogeny from a coalescent genealogy
#' @param phy genealogy
#' @param mu mu
#' @param omega omega
#' @return a mutation scaled phylogeny
simulate_mut_arc <- function(phy, mu, omega)
{
    phy_mut <- phy

    theta <- mu / omega
    p <- 1/(1+omega)

    sim_mut <- function(x) if(x > 1e-8) rnbinom(1, theta * x, p) else 0 
    phy_mut$edge.length <- sapply(phy_mut$edge.length, sim_mut) 
    
    #phy_mut$edge.length <- rnbinom(length(phy_mut$edge.length), theta * phy_mut$edge.length, rep(p,length(phy_mut$edge.length)))
    
    #Randomise polytomy splits
    phy_mut <- di2multi(phy_mut)
    phy_mut <- multi2di(phy_mut,random=TRUE,equiprob=FALSE)

    return(phy_mut)
}

validate_input <- function(samp_times, n_samp)
{
    stopifnot("sampling times must be in ascending order"=!is.unsorted(samp_times))
    stopifnot("samp_times must have same length as n_samp"=length(samp_times) == length(n_samp))
    stopifnot("n_samp must be integer valued"=is.integer(n_samp))
    stopifnot("n_samp must be positive"=all(n_samp > 0))
}