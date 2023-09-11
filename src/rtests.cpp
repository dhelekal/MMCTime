#include <Rcpp.h>
#include <algorithm>
#include <tuple>
#include <vector>
#include "state_manip.h"
#include "tree_data.h"
#include "par_data.h"
#include "lambda.h"
#include "probs.h"


/*
 This file contains some extra functions to simplify proposal ratio testing for non-RWM moves
*/

// [[Rcpp::export]]
double taus_to_times_logJ(Rcpp::List phydata, Rcpp::List mod_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);   
    return mmctime::tau_log_jacobian(bt, st); 
}

// [[Rcpp::export]]
double times_to_taus_logJ(Rcpp::List phydata, Rcpp::List mod_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);   
    return mmctime::tau_inv_log_jacobian(bt, st); 
}

// [[Rcpp::export]]
double kingman_lp(Rcpp::List phydata, Rcpp::List mod_state, double nu)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);  

    const mmctime::kingman coal{nu};

    auto events = mmctime::convert_state(st,bt);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
    return mmctime::lambda_lp(events, bt.n_tip(), coal);
}

// [[Rcpp::export]]
double km_beta_lp(Rcpp::List phydata, Rcpp::List mod_state, double nu, double phi, double alpha)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata); 

    const mmctime::km_beta_lambda coal{phi, nu, 1.0+alpha};

    auto events = mmctime::convert_state(st,bt);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
    return mmctime::lambda_lp(events, bt.n_tip(), coal);
}

// [[Rcpp::export]]
double ds_lp(Rcpp::List phydata, Rcpp::List mod_state, double nu, double phi)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata); 

    const mmctime::ds_lambda coal{phi, nu};

    auto events = mmctime::convert_state(st,bt);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
    return mmctime::lambda_lp(events, bt.n_tip(), coal);     
}

// [[Rcpp::export]]
double beta_lp(Rcpp::List phydata, Rcpp::List mod_state, double nu, double alpha)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);  

    const mmctime::beta_lambda coal{nu, 1.0+alpha};

    auto events = mmctime::convert_state(st,bt);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
    return mmctime::lambda_lp(events, bt.n_tip(), coal);
}