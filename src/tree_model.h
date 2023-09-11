#ifndef TREE_MODEL
#define TREE_MODEL

#include <Rcpp.h>
#include <stdexcept>
#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"
#include <iostream>
#include <vector>


namespace mmctime
{
    template<typename T, typename U>
    struct tree_model
    {
        tree_model(par_state ps, T coal, U obs) :
            pars(ps), coal_impl(coal), obs_impl(obs) {};

        double transform_pars();
        double lp_coal(const btree& btree, const mc_state& mc_state) const;
        double lp_obs(const btree& btree, const mc_state& mc_state, const branch_data_IS& obs) const;
        std::vector<double> summaries(const btree& btree, const mc_state& mc_state) const;

        int n_pars_coal() const;
        int n_pars_obs() const;
        par_state pars;
        T coal_impl;
        U obs_impl;
    };  
} 

template<typename T, typename U>
double mmctime::tree_model<T,U>::transform_pars()
{
    const int n = obs_impl.n_pars;
    const int m = coal_impl.n_pars;

    if (!(pars.size() == n + m)) throw 
        std::invalid_argument("Parameter vector size doesn't match!");

    return obs_impl.transform_pars(pars.create_view(0, n)) + coal_impl.transform_pars(pars.create_view(n, m));
}

template<typename T, typename U>
double mmctime::tree_model<T,U>::lp_coal(const btree& btree, const mc_state& mc_state) const
{
    const int n = obs_impl.n_pars;
    const int m = coal_impl.n_pars;
    if (!(pars.size() == n + m)) throw 
        std::invalid_argument("Parameter vector size doesn't match!");
    const double lp_coal = coal_impl.prior_lp_tree(btree, mc_state, pars.create_view(n, m));
    return lp_coal;
}

template<typename T, typename U>
double mmctime::tree_model<T,U>::lp_obs(const btree& btree, const mc_state& mc_state, const branch_data_IS& obs) const
{
    const int n = obs_impl.n_pars;
    const int m = coal_impl.n_pars;
    if (!(pars.size() == n+m)) throw 
        std::invalid_argument("Parameter vector size doesn't match!");
    
    return obs_impl.obs_lp_branch(btree, mc_state, obs, pars.create_view(0, n));
}

template<typename T, typename U>
std::vector<double> mmctime::tree_model<T,U>::summaries(const btree& btree, const mc_state& mc_state) const
{
    const int n = obs_impl.n_pars;
    const int m = coal_impl.n_pars;
    if (!(pars.size() == n + m)) throw 
        std::invalid_argument("Parameter vector size doesn't match!");
    
    return coal_impl.summaries(btree, mc_state, pars.create_view(n, m));
}

template<typename T, typename U>
int mmctime::tree_model<T,U>::n_pars_coal() const
{
    return coal_impl.n_pars;
}

template<typename T, typename U>
int mmctime::tree_model<T,U>::n_pars_obs() const
{
    return obs_impl.n_pars;
}

#endif