#include "rwm_move.h"
#include "state_manip.h"
#include "prob_utils.h"
#include <Rcpp.h>
#include <cmath>

double mmctime::tau_prop::rwm_prop(const btree& btree, const mc_state& curr_state, mc_state& next_state) const
{
    const int root_idx = btree.root_idx();
    next_state.copy_from(curr_state);
    int j = 0;
    for (int i = btree.n_tip(); i < btree.n_node(); i++)
    {
        if (i == root_idx)
        {
            const double eps = R::rnorm(0.0, 1.0) * sd.at(j);
            const double tau_prev = next_state.tau_at(i);
            next_state.tau_at(i) = tau_prev + eps;
        } else if (curr_state.q_at(i) == 0) 
        {
            const double eps = R::rnorm(0.0, 1.0) * sd.at(j);
            const double tau_prev = next_state.tau_at(i);
            next_state.tau_at(i) = std::fabs(tau_prev + eps);
        } 
        j++;
    }

    times_from_taus(btree, next_state, btree.root_idx());
    return 0.0;
}

double mmctime::par_prop::rwm_prop(const par_state& curr_state, par_state& next_state) const
{
    next_state.copy_from(curr_state);
    auto& next = next_state.unconstr(); 

    for (int i = 0; i < next_state.size(); i++)
    {
        next.at(i) += R::rnorm(0.0,sd.at(i));
    }
   
    return 0.0;
}
