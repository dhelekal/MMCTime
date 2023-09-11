#ifndef OBS_MODELS
#define OBS_MODELS

#include <Rcpp.h>
#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"

namespace mmctime
{
    struct arc_model
    {
        double transform_pars(par_state ps);
        double obs_lp_branch(const btree& btree, const mc_state& mc_state, const branch_data_IS& obs, const par_state ps) const;
        int n_pars = 2;
    };
}
#endif