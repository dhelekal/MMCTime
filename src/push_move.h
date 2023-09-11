#ifndef PUSH_MOVE
#define PUSH_MOVE

#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"
#include <tuple>

namespace mmctime
{   
    class mm_prop2
    {
    public: 
        double prop(const btree& btree, const mc_state& curr_state, mc_state& next_state) const;
        
        double move_bp(const btree& btree, const mc_state& prev_state, mc_state& next_state, int node_idx) const;
        double qr_bp(const btree& btree, const mc_state& prev_state, const mc_state& next_state, int node_idx) const;
        double lp_node(const btree& btree, const mc_state& state, int node_idx) const;
        
        double _sd;
        branch_data_IS _br_data;

    private:

        std::tuple<double, int> sample_tau_idx(const mc_state& mcs, const btree& btree) const;

        double tau_lwfunc(const mc_state& mcs, const btree& btree, int n_idx) const;
        double tau_lwsum(const mc_state& mcs, const btree& btree) const;
        
        bool is_mm_adj(const btree& btree,
            const mc_state& new_state,
            int node_idx) const;

        double merge_up_rec(const btree& btree,
            mc_state& new_state,
            int node_idx) const;

        double split_down_rec(const btree& btree,
            mc_state& new_state,
            bool first_node,
            int node_idx) const;

        double expand_one(const btree& btree,
            mc_state& new_state,
            int node_idx) const;
        
        double collapse_one(const btree& btree,
            mc_state& new_state,
            int node_idx) const;

        double sample_tau() const;
        double prop_lp_tau(double tau) const;
    };
}

#endif