//probs.cpp
#include "probs.h"
#include "prob_utils.h"
#include "tree_data.h"
#include "par_data.h"
#include <tuple>
#include <functional>
#include <Rcpp.h>

double mmctime::tau_log_jacobian(const btree& btree, const mc_state& mc_state)
{

    double log_jac = 0.0;
    for (int i = btree.n_tip(); i < btree.n_node(); i++)
    {
        if (i == btree.root_idx())
        {
            const double tau = mc_state.tau_at(i);
            log_jac += tau;
        } else if (mc_state.q_at(i) == 0)
        {
            const double pa_t = mc_state.t_at(btree.get_parent(i));
            const double sb = btree.get_sbound(i);
            const double tau = mc_state.tau_at(i);
            log_jac += std::log(pa_t - sb) - tau;
        }
    }
    return log_jac;
}


double mmctime::tau_inv_log_jacobian(const btree& btree, const mc_state& mc_state)
{
    double log_jac = 0.0;
    for (int i = btree.n_tip(); i < btree.n_node(); i++)
    {
        if (i == btree.root_idx())
        {
            const double t = mc_state.t_at(i);
            const double sb = btree.get_sbound(i);
            log_jac += -std::log(t - sb);
        } else if (mc_state.q_at(i) == 0) 
        {
            const double sb = btree.get_sbound(i);
            const double t = mc_state.t_at(i);
            log_jac += -std::log(t - sb);
        }
    }
    return log_jac;
}