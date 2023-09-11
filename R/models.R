tree_model <- function(model)
{
    out <- NA
    if (model == "km_beta")
    {
        out <- list(transform_pars = transform_kmb_arc,
            obs_lp = prob_kmb_arc, 
            coal_lp = coal_lp_kmb_arc,
            summaries = summaries_kmb_arc,
            n_par_tree = 3,
            n_par_obs = 2,
            names_par_tree = c("nu", "phi", "alpha"),
            names_par_obs = c("mu", "omega"),
            names_summaries = c("mm_count", "max_msize")
        )
    } else if (model == "beta")
    {
        out <- list(transform_pars = transform_beta_arc,
            obs_lp = prob_beta_arc, 
            coal_lp = coal_lp_beta_arc,
            summaries = summaries_beta_arc,
            n_par_tree = 2,
            n_par_obs = 2,
            names_par_tree = c("nu", "alpha"),
            names_par_obs = c("mu", "omega"),
            names_summaries = c("mm_count","max_msize")
        )
    } else if (model =="durret-schweinsberg")
    {
        out <- list(transform_pars = transform_ds_arc,
            obs_lp = prob_ds_arc, 
            coal_lp = coal_lp_ds_arc,
            summaries = summaries_ds_arc,
            n_par_tree = 2,
            n_par_obs = 2,
            names_par_tree = c("nu", "phi"),
            names_par_obs = c("mu", "omega"),
            names_summaries = c("mm_count","max_msize")
        )
    } else if (model == "kingman") 
    {
        out <- list(transform_pars = transform_km_arc,
            obs_lp = prob_km_arc, 
            coal_lp = coal_lp_km_arc,
            summaries = summaries_km_arc,
            n_par_tree = 1,
            n_par_obs = 2,
            names_par_tree = c("nu"),
            names_par_obs = c("mu", "omega"),
            names_summaries = c()
        )
    } else
    {
        stop("Unrecognised Model.")
    }
    return(out)
}

MODEL_OPTS <- function() return(c("km_beta", "beta", "durret-schweinsberg", "kingman"))