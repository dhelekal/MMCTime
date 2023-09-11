#include <Rcpp.h>
#include "sim_lambda.h"
#include "tree_init.h"
#include "state_manip.h"
#include "tree_data.h"
#include "branch_data.h"
#include "par_data.h"
#include "lambda.h"
#include <tuple>
#include <vector>

// [[Rcpp::export]]
Rcpp::List sim_km_beta(Rcpp::NumericVector samp_times, Rcpp::IntegerVector n_samp, double phi, double nu, double alpha)
{   
    mmctime::km_beta_lambda coal{phi, nu, 1.0+alpha};
    return  mmctime::sim_lambda_direct(samp_times, n_samp, coal);
}

// [[Rcpp::export]]
Rcpp::List sim_beta(Rcpp::NumericVector samp_times, Rcpp::IntegerVector n_samp, double nu, double alpha)
{
    mmctime::beta_lambda coal{nu,1.0+alpha};
    return  mmctime::sim_lambda_direct(samp_times, n_samp, coal);
}

// [[Rcpp::export]]
Rcpp::List sim_kingman(Rcpp::NumericVector samp_times, Rcpp::IntegerVector n_samp, double nu)
{
    mmctime::kingman coal{nu};
    return  mmctime::sim_lambda_direct(samp_times, n_samp, coal);
}

// [[Rcpp::export]]
Rcpp::List sim_ds(Rcpp::NumericVector samp_times, Rcpp::IntegerVector n_samp, double phi, double nu)
{   
    mmctime::ds_lambda coal{phi, nu};
    return mmctime::sim_lambda_direct(samp_times, n_samp, coal);
}

// [[Rcpp::export]]
void initialise_times(Rcpp::List phydata, Rcpp::List mod_state, Rcpp::IntegerVector branch_data, double mu)
{
    const int n_node = phydata["n_node"];
    const int n_tip = phydata["n_tip"];
    const int root_br_pos = phydata["root_br_pos"];

    const Rcpp::IntegerMatrix& topo_mat = phydata["topo_mat"];
    Rcpp::NumericVector node_times = mod_state["times"];

    const int root_idx = n_node-1;
    
    mmctime::init_times(topo_mat, branch_data, n_node, n_tip, mu, root_br_pos, root_idx, node_times);
}

// [[Rcpp::export]]
Rcpp::List duplicate_mc_state(Rcpp::List mod_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    return st.copy().as_list();
}

// [[Rcpp::export]]
void mc_state_copy_to(Rcpp::List src, Rcpp::List dest)
{
    const auto st_src = mmctime::mc_state::from_list(src);
    auto st_dest = mmctime::mc_state::from_list(dest);
    st_dest.copy_from(st_src);
}

// [[Rcpp::export]]
Rcpp::List duplicate_par_state(Rcpp::List mod_state)
{
    const auto ps = mmctime::par_state::from_list(mod_state);
    return ps.copy().as_list();
}

// [[Rcpp::export]]
void par_state_copy_to(Rcpp::List src, Rcpp::List dest)
{
    const auto ps_src = mmctime::par_state::from_list(src);
    auto ps_dest = mmctime::par_state::from_list(dest);
    ps_dest.copy_from(ps_src);
}

// [[Rcpp::export]]
Rcpp::List duplicate_tree(Rcpp::List phydata)
{
    const auto tr = mmctime::btree::from_list(phydata);
    return tr.copy().as_list();
}

// [[Rcpp::export]]
void tree_copy_to(Rcpp::List src, Rcpp::List dest)
{
    const auto tr_src = mmctime::btree::from_list(src);
    auto tr_dest = mmctime::btree::from_list(dest);
    tr_dest.copy_from(tr_src);
}

// [[Rcpp::export]]
void times_from_taus(Rcpp::List phydata, Rcpp::List mod_state)
{
    auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);
    
    mmctime::times_from_taus(bt, st, bt.root_idx());
}

//' @export
// [[Rcpp::export]]
void taus_from_times(Rcpp::List phydata, Rcpp::List mod_state)
{
    auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);

    mmctime::taus_from_times(bt, st, bt.root_idx());
}

// [[Rcpp::export]]
bool validate_times(Rcpp::List phydata, Rcpp::List mod_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);

    return mmctime::validate_times(bt, st, bt.root_idx());
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix build_topo_mat(Rcpp::IntegerMatrix br_mat, int n_node, int root_pos, double cutoff)
{
    return mmctime::build_topo_mat(br_mat, n_node, root_pos, cutoff);
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix build_branch_data(Rcpp::IntegerMatrix br_mat, int n_node, int root_pos)
{
    return mmctime::build_branch_mat(br_mat, n_node, root_pos);
}

// [[Rcpp::export]]
Rcpp::NumericVector find_sbounds(Rcpp::IntegerMatrix br_mat, Rcpp::NumericVector tip_times, int n_tip, int root_idx)
{
    return mmctime::create_sbounds(br_mat, tip_times, n_tip, root_idx);
}

//' @export
// [[Rcpp::export]]
void binary_to_mm(Rcpp::List phydata, Rcpp::List mod_state, double tol=1e-8)
{
    auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata);
    mmctime::squash_brs(bt, tol, st, bt.root_idx());
    mmctime::times_from_taus(bt, st, bt.root_idx());
}

// [[Rcpp::export]]
Rcpp::IntegerVector find_pivots(Rcpp::List phydata, Rcpp::IntegerVector branch_data)
{
    std::vector<int> pivots;
    const auto bt = mmctime::btree::from_list(phydata);
    const mmctime::branch_data_IS dat(branch_data);

    bool root_fixed = bt.root_fixed() > 0;

    const int r_idx = bt.root_idx();
    int r_c1, r_c2;
    std::tie(r_c1, r_c2) = bt.get_children(r_idx);

    //If root position is fixed enabled polytomy interchanges through the root
    //Disabled otherwise to keep proposal ratios simple as RootNNI shadows the move

    for (int i = bt.n_tip(); i < bt.n_node(); i++)
    {        
        if ((root_fixed || (i != r_c1 && i != r_c2)) && i != r_idx && dat.is_poly(bt.pa_edge(i)))
        {
            pivots.push_back(i);
        }
    }
    return Rcpp::wrap(pivots);
}

// [[Rcpp::export]]
Rcpp::List as_lambda_events(Rcpp::List phydata, Rcpp::List mod_state)
{
    const auto st = mmctime::mc_state::from_list(mod_state);
    const auto bt = mmctime::btree::from_list(phydata); 

    const auto evts = mmctime::convert_state(st,bt);
    const int sz = evts.size();


    Rcpp::NumericVector evt_times(sz);
    Rcpp::IntegerVector delta_At(sz);

    for (int i = 0; i < sz; i++)
    {
        auto d = evts.at(i);
        evt_times.at(i) = std::get<0>(d);
        delta_At.at(i) = std::get<1>(d); 
    }

    return Rcpp::List::create(Rcpp::Named("event_times") = evt_times,
        Rcpp::Named("delta_At") = delta_At);
}
