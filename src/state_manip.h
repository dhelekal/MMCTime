#ifndef STATE_MANIP
#define STATE_MANIP

#include "tree_data.h"
#include "par_data.h"
#include <vector>
#include <tuple>

namespace mmctime
{
    std::vector<std::tuple<double,int>> convert_state(const mc_state& mc_state, const btree& btree);

    void taus_from_times(const btree& btree, mc_state& mc_state, int node_idx);
    void times_from_taus(const btree& btree, mc_state& mc_state, int node_idx);

    int q_count(const mc_state& mc_state);

    double mm_lower_bound(const mc_state& mc_state,
        const btree& btree,
        int node_idx);

    int mm_count(const mc_state& mc_state,
        const btree& btree,
        int node_idx);

    int q_rec(const mc_state& mc_state,
        const btree& btree,
        int node_idx);

    void mm_update_taus(const btree& btree,
        mc_state& mc_state,
        int node_idx,
        double t);

    void squash_brs(const btree& btree,
        double tol,
        mc_state& mc_state,
        int node_idx);

    bool validate_times(const btree& btree,
        const mc_state& mc_state,
        int node_idx);

    void pars_from_unconstrained(par_state& par_state);
    void unconstrained_from_pars(par_state& par_state);
}


#endif