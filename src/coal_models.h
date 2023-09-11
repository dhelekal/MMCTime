#ifndef COAL_MODELS
#define COAL_MODELS

#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"
#include <vector>

namespace mmctime
{   
    struct km_beta_model
    {
        double transform_pars(par_state ps) const;
        double prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        std::vector<double> summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        int n_pars = 3;
    };

    struct ds_model
    {
        double transform_pars(par_state ps) const;
        double prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        std::vector<double> summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        int n_pars = 2;
    };

    struct beta_model
    {
        double transform_pars(par_state ps) const;
        double prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        std::vector<double> summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        int n_pars = 2;
    };

    struct km_model
    {
        double transform_pars(par_state ps) const;
        double prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        std::vector<double> summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const;
        int n_pars = 1;
    }; 
}
#endif