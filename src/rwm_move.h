#ifndef RWM_MOVE
#define RWM_MOVE

#include "tree_data.h"
#include "par_data.h"
#include <Rcpp.h>

namespace mmctime
{
    struct tau_prop
    {   
        Rcpp::NumericVector sd;
        double rwm_prop(const btree& btree, const mc_state& curr_state, mc_state& next_state) const;
    };

    struct par_prop
    {
        Rcpp::NumericVector sd; 
        double rwm_prop(const par_state& curr_state, par_state& next_state) const;
    };
}

#endif