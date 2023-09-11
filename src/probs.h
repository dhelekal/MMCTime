//probs.h
#ifndef PROBS_H // include guard
#define PROBS_H

#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"
#include "state_manip.h"
#include "prob_utils.h"
#include "lambda.h"
#include <Rcpp.h>
#include <tuple>
#include <algorithm>
namespace mmctime 
{
    double tau_log_jacobian(const btree& btree, const mc_state& mc_state);
    double tau_inv_log_jacobian(const btree& btree, const mc_state& mc_state);

    double prob_mm_topo(const btree& btree, const mc_state& mc_state);

    template <typename T>
    double prior_lp_lambda(const btree& btree, const mc_state& mc_state, T lambda_coal);

    template <typename T>
    double prior_lp_mix(const btree& btree, const mc_state& mc_state, double mix_prob, T mix_coal, kingman km_coal);

    template <typename T>
    double coal_mix_ratio(const btree& btree, const mc_state& mc_state,  double mix_prob, T mix_coal, kingman km_coal);

    template <typename T>
    double obs_lp(const btree& btree, const mc_state& mc_state, T obs);
}

template <typename T>
double mmctime::prior_lp_lambda(const btree& btree, const mc_state& mc_state, T lambda_coal)
{
    double out = 0.0; 
    
    auto events = convert_state(mc_state, btree);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});

    const int n_tip = btree.n_tip();

    out += lambda_lp(events, n_tip, lambda_coal);
    out += orthant_correction(events);
    out += tau_log_jacobian(btree, mc_state);
    return out;
}

template <typename T>
double mmctime::prior_lp_mix(const btree& btree, const mc_state& mc_state, double mix_prob, T mix_coal, kingman km_coal)
{
    double out = 0.0; 
    const double a = std::log(mix_prob);
    auto events = convert_state(mc_state, btree);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
    
    const int max_msize = 1-std::get<1>(*std::max_element(events.begin(),
        events.end(), 
        [] (auto a, auto b) {return std::get<1>(a) > std::get<0>(b);}));

    const int n_tip = btree.n_tip();

    double lp_km = a;
    if (max_msize > 2)
    {
        lp_km += -INFINITY;
    } else 
    {   
        lp_km += lambda_lp(events, n_tip, km_coal);
    } 
    const double lp_mix = log_1m_exp(a) + lambda_lp(events, n_tip, mix_coal);
    out += log_sum_exp(lp_km, lp_mix);
    out += tau_log_jacobian(btree, mc_state);
    out += orthant_correction(events);

    return out;
}

template <typename T>
double mmctime::coal_mix_ratio(const btree& btree, const mc_state& mc_state, double mix_prob, T mix_coal, kingman km_coal)
{
    const double a = std::log(mix_prob);
    
    auto events = convert_state(mc_state, btree);
    std::sort(events.begin(), events.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
    
    const int max_msize = 1 - std::get<1>(*std::max_element(events.begin(),
        events.end(), 
        [] (auto a, auto b) {return std::get<1>(a) > std::get<0>(b);}));

    const int n_tip = btree.n_tip();

    double lp_km = a;
    if (max_msize > 2)
    {
        lp_km += -INFINITY;
    } else 
    {   
        lp_km += lambda_lp(events, n_tip, km_coal);
    } 

    const double lp_mix = log_1m_exp(a) + lambda_lp(events, n_tip, mix_coal);
    const double lp_tot = log_sum_exp(lp_km, lp_mix);

    return lp_mix-lp_tot;
}

template <typename T>
double mmctime::obs_lp(const btree& btree, const mc_state& mc_state, const branch_data_IS& data, T obs)
{
    const int r_n = btree.root_idx();
    int r_c1, r_c2;
    std::tie(r_c1, r_c2) = btree.get_children(r_n);

    double out = 0.0;

    //Non root adjacent edge probs
    for (int i = 0; i < btree.n_node(); i++)
    {
        if (i == r_n || i == r_c1 || i == r_c2) continue;
        
        const int br_above = btree.pa_edge(i);

        const int n_mut = data.get_branch_muts(br_above);
        if (n_mut > 0 && mc_state.q_at(i) == 1)
        {
            out += -INFINITY;
        } else if (mc_state.q_at(i) == 0)
        {
            const int pa_n = btree.get_parent(i);
            const double br_len = mc_state.t_at(pa_n) - mc_state.t_at(i);
            out += obs.lp(n_mut, br_len);
        } 
    }

    double lp_root = 0;
    if (btree.root_fixed())
    {
        // If root branch is fixed then the tree is rooted a-priori
        int rcs[2] {r_c1, r_c2};
        for (int i : rcs)
        {
            const double brlen = mc_state.t_at(r_n) - mc_state.t_at(i);
            const int br_above = btree.pa_edge(i);
            const int n_mut = data.get_branch_muts(br_above);
            lp_root += obs.lp(n_mut, brlen);
        }
    } else
    {
        // If root branch is not fixed marginalise out position of root on root branch
        const double r_brlen1 = mc_state.t_at(r_n) - mc_state.t_at(r_c1);
        const double r_brlen2 = mc_state.t_at(r_n) - mc_state.t_at(r_c2);
        const int root_mut = data.get_branch_muts(btree.root_br_pos());
        lp_root += obs.lp(root_mut, r_brlen1 + r_brlen2);
    }

    out += lp_root;
    return out;
}


#endif

