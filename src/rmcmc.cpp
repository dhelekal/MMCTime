#include <Rcpp.h>
#include "tree_data.h"
#include "par_data.h"
#include "rwm_move.h"
#include "push_move.h"
#include "root_move.h"
#include "topo_move.h"

// [[Rcpp::export]]
double rwm_tau_move(Rcpp::List phydata, Rcpp::List curr_mc_state, Rcpp::List prop_mc_state, Rcpp::NumericVector sd)
{
    const auto bt = mmctime::btree::from_list(phydata);
    const auto st_curr = mmctime::mc_state::from_list(curr_mc_state);
    auto st_prop = mmctime::mc_state::from_list(prop_mc_state);

    mmctime::tau_prop prop{sd};

    return prop.rwm_prop(bt, st_curr, st_prop);
}

// [[Rcpp::export]]
double rwm_par_move(Rcpp::List curr_par_state, Rcpp::List prop_par_state, Rcpp::NumericVector sd)
{
    const auto ps_curr = mmctime::par_state::from_list(curr_par_state);
    auto ps_prop = mmctime::par_state::from_list(prop_par_state);

    mmctime::par_prop prop{sd};

    return prop.rwm_prop(ps_curr, ps_prop);
}

// [[Rcpp::export]]
double push_mm_move(Rcpp::List phydata, Rcpp::List curr_state, Rcpp::List prop_state, double sd, Rcpp::IntegerVector branch_data)
{
    const auto bt = mmctime::btree::from_list(phydata);
    const auto st_curr = mmctime::mc_state::from_list(curr_state);
    
    mmctime::branch_data_IS br_data(branch_data);

    auto st_prop = mmctime::mc_state::from_list(prop_state);

    mmctime::mm_prop2 prop{sd, br_data};

    return prop.prop(bt, st_curr, st_prop);
}

// [[Rcpp::export]]
double topo_move(Rcpp::List curr_phydata, Rcpp::List curr_state, Rcpp::List prop_phydata, Rcpp::List prop_state, Rcpp::IntegerVector branch_data, Rcpp::IntegerVector pivots)
{
    const auto bt_curr = mmctime::btree::from_list(curr_phydata);
    const auto st_curr = mmctime::mc_state::from_list(curr_state);
    
    auto bt_prop = mmctime::btree::from_list(prop_phydata);
    auto st_prop = mmctime::mc_state::from_list(prop_state);

    const mmctime::branch_data_IS br_data(branch_data);

    mmctime::topo_prop prop{pivots, br_data};

    return prop.prop(bt_curr, st_curr, bt_prop, st_prop);
}

// [[Rcpp::export]]
double root_move(Rcpp::List curr_phydata, Rcpp::List curr_state, Rcpp::List prop_phydata, Rcpp::List prop_state)
{
    const auto bt_curr = mmctime::btree::from_list(curr_phydata);
    const auto st_curr = mmctime::mc_state::from_list(curr_state);
    
    auto bt_prop = mmctime::btree::from_list(prop_phydata);
    auto st_prop = mmctime::mc_state::from_list(prop_state);

    mmctime::root_prop prop;

    return prop.prop(bt_curr, st_curr, bt_prop, st_prop);
}


