#include <Rcpp.h>
#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"
#include "tree_model.h"
#include "coal_models.h"
#include "obs_models.h"

typedef mmctime::km_beta_model kmb;
typedef mmctime::ds_model ds;
typedef mmctime::beta_model beta;
typedef mmctime::km_model km;
typedef mmctime::arc_model arc;

template<typename MOD_COAL, typename MOD_OBS>
inline double comp_lp_tmp(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state, Rcpp::IntegerVector branch_data)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata); 
    const auto ps = mmctime::par_state::from_list(par_state);
    const mmctime::branch_data_IS br_data(branch_data);

    MOD_COAL mod_coal;
    MOD_OBS mod_obs;

    mmctime::tree_model<MOD_COAL, MOD_OBS> tm(ps, mod_coal, mod_obs);
    return tm.lp_obs(bt, st, br_data);
}

template<typename MOD_COAL, typename MOD_OBS>
inline double comp_coal_lp_tmp(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata); 
    const auto ps = mmctime::par_state::from_list(par_state);

    MOD_COAL mod_coal;
    MOD_OBS mod_obs;

    mmctime::tree_model<MOD_COAL, MOD_OBS> tm(ps, mod_coal, mod_obs);
    return tm.lp_coal(bt, st);
}

template<typename MOD_COAL, typename MOD_OBS>
inline double transf_pars_tmp(Rcpp::List par_state)
{
    auto ps = mmctime::par_state::from_list(par_state);

    MOD_COAL mod_coal;
    MOD_OBS mod_obs;
    
    mmctime::tree_model<MOD_COAL, MOD_OBS> tm(ps, mod_coal, mod_obs);
    return tm.transform_pars();
}

template<typename MOD_COAL, typename MOD_OBS>
inline Rcpp::NumericVector comp_summaries_tmp(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata); 
    const auto ps = mmctime::par_state::from_list(par_state);

    MOD_COAL mod_coal;
    MOD_OBS mod_obs;

    mmctime::tree_model<MOD_COAL, MOD_OBS> tm(ps, mod_coal, mod_obs);
    return Rcpp::wrap(tm.summaries(bt, st));
}

/*
BEGIN KM-BETA MODEL FUNCS
*/
// [[Rcpp::export]]
double transform_kmb_arc(Rcpp::List par_state)
    {return transf_pars_tmp<kmb, arc>(par_state);}

// [[Rcpp::export]]
double prob_kmb_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state, Rcpp::IntegerVector branch_data)
    {return comp_lp_tmp<kmb, arc>(phydata, mod_state, par_state, branch_data);}

// [[Rcpp::export]]
double coal_lp_kmb_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_coal_lp_tmp<kmb, arc>(phydata, mod_state, par_state);}

// [[Rcpp::export]]
Rcpp::NumericVector summaries_kmb_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_summaries_tmp<kmb, arc>(phydata, mod_state, par_state);}
/*
END KM-BETA MODEL FUNCS
*/

/*
BEGIN DS MODEL FUNCS
*/
// [[Rcpp::export]]
double transform_ds_arc(Rcpp::List par_state)
    {return transf_pars_tmp<ds, arc>(par_state);}

// [[Rcpp::export]]
double prob_ds_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state, Rcpp::IntegerVector branch_data)
    {return comp_lp_tmp<ds, arc>(phydata, mod_state, par_state, branch_data);}

// [[Rcpp::export]]
double coal_lp_ds_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_coal_lp_tmp<ds, arc>(phydata, mod_state, par_state);}
// [[Rcpp::export]]
Rcpp::NumericVector summaries_ds_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_summaries_tmp<ds, arc>(phydata, mod_state, par_state);}
/*
END DS MODEL FUNCS
*/

/*
BEGIN BETA MODEL FUNCS
*/
// [[Rcpp::export]]
double transform_beta_arc(Rcpp::List par_state)
    {return transf_pars_tmp<beta, arc>(par_state);}

// [[Rcpp::export]]
double prob_beta_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state, Rcpp::IntegerVector branch_data)
    {return comp_lp_tmp<beta, arc>(phydata, mod_state, par_state, branch_data);}

// [[Rcpp::export]]
double coal_lp_beta_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_coal_lp_tmp<beta, arc>(phydata, mod_state, par_state);}
// [[Rcpp::export]]
Rcpp::NumericVector summaries_beta_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_summaries_tmp<beta, arc>(phydata, mod_state, par_state);}
/*
END BETA MODEL FUNCS
*/

/*
BEGIN KM MODEL FUNCS
*/
// [[Rcpp::export]]
double transform_km_arc(Rcpp::List par_state)
    {return transf_pars_tmp<km, arc>(par_state);}

// [[Rcpp::export]]
double prob_km_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state, Rcpp::IntegerVector branch_data)
    {return comp_lp_tmp<km, arc>(phydata, mod_state, par_state, branch_data);}

// [[Rcpp::export]]
double coal_lp_km_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_coal_lp_tmp<km, arc>(phydata, mod_state, par_state);}
// [[Rcpp::export]]
Rcpp::NumericVector summaries_km_arc(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::List par_state)
    {return comp_summaries_tmp<km, arc>(phydata, mod_state, par_state);}
/*
END KM MODEL FUNCS
*/