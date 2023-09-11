alpha_prior <- function(mod, alpha_shape1=get_default("alpha_shape1", mod), alpha_shape2=get_default("alpha_shape2", mod), ...) {return(function(y) dbeta(y, alpha_shape1, alpha_shape2, log=T))}
nu_prior <- function(mod, nu_logmean=get_default("nu_logmean", mod), nu_logsd=get_default("nu_logsd", mod), ...) {return(function(y) dlnorm(y, nu_logmean, nu_logsd, log=T))}
phi_prior <- function(mod, phi_shape1=get_default("phi_shape1", mod), phi_shape2=get_default("phi_shape2", mod), ...) {return(function(y) dbeta(y, phi_shape1, phi_shape2, log=T))}
mu_prior <- function(mod, mu_shape=get_default("mu_shape", mod), mu_scale=get_default("mu_scale", mod), ...) {return(function(y) dgamma(y, shape=mu_shape, scale=mu_scale, log=T))}
omega_prior <- function(mod, omega_mu=get_default("omega_mu", mod), omega_sd=get_default("omega_sd", mod), ...) 
{
    a <- log(2.0)
    return(function(y) a + dnorm(y, omega_mu, omega_sd, log=T))
}

prior_builder <- function(NAME, mod, ...)
{
    lookup <- list(nu=nu_prior, phi=phi_prior, alpha=alpha_prior, mu=mu_prior, omega=omega_prior)
    return(lookup[[NAME]](mod, ...))
}

get_default <- function(pn, mod)
{
    if (pn == "alpha_shape1")
    {
        return(switch(mod, km_beta=PPAR_DEFAULTS[["alpha_shape1_km_beta"]],
            beta=PPAR_DEFAULTS[["alpha_shape1_beta"]],
            stop("Unrecognised model / parameter")))
    }
    else if (pn == "alpha_shape2")
    {
        return(switch(mod, km_beta=PPAR_DEFAULTS[["alpha_shape2_km_beta"]],
            beta=PPAR_DEFAULTS[["alpha_shape2_beta"]],
            stop("Unrecognised model / parameter")))
    }
    else
    {
        return(PPAR_DEFAULTS[[pn]])
    }
}

#' Prior hyperparameter names
#' @export 
PPAR_NAMES <- c(
    "alpha_shape1", 
    "alpha_shape2", 
    "nu_logmean", 
    "nu_logsd", 
    "phi_shape1",
    "phi_shape2", 
    "mu_shape", 
    "mu_scale", 
    "omega_mu", 
    "omega_sd")

#' Prior hyperparameter default values
#' @export
PPAR_DEFAULTS <- list(   
    alpha_shape1_km_beta = 1,
    alpha_shape2_km_beta = 2,
    alpha_shape1_beta = 3,
    alpha_shape2_beta = 1,
    nu_logmean=0, 
    nu_logsd=4, 
    phi_shape1=1, 
    phi_shape2=3,
    mu_shape=2, 
    mu_scale=8, 
    omega_mu=0, 
    omega_sd=2
)

make_par_prior <- function(PAR_NAMES, mod, ...)
{
    la <- list(...)
    ln <- names(la)
    for (pn in ln)
    {
        if (!(pn %in% PPAR_NAMES))
            warning(paste0("Unrecognised prior hyperparameter: (", pn,"=",la[[pn]], ")"))
    }

    priors <- lapply(PAR_NAMES, function(x) prior_builder(x, mod, ...))
    names(priors) <- PAR_NAMES

    n_par <- length(PAR_NAMES)

    compute_prior_lp <- function(x)
    {
        lp <- 0 
        for (i in 1:n_par)
        {
            lp <- lp + priors[[PAR_NAMES[[i]]]](x[i])
        }
        return(lp)
    }
    return(list(par_priors=priors, prior_lp=compute_prior_lp))
}