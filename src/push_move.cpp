#include "push_move.h"
#include "tree_data.h"
#include "state_manip.h"
#include "prob_utils.h"
#include <functional>
#include <stdexcept>
#include <Rcpp.h>

double mmctime::mm_prop2::prop(const btree& btree, const mc_state& curr_state, mc_state& next_state) const
{
    next_state.copy_from(curr_state);
    double out = 0.0;

    double n_lp; int n_idx;
    std::tie(n_lp, n_idx) = sample_tau_idx(next_state, btree);
    
    // Lp of proposal selecting node
    out += -n_lp;

    const bool is_merged = (curr_state.q_at(n_idx) == 1);
    const bool below_mm = is_mm_adj(btree, curr_state, n_idx);

    const double loghalf = std::log(0.5);

    if (!below_mm && !is_merged)
    {
        // Node is not merged as a part of an MM & is not directly above a multi merger
        // Eligible for branching process merge move
        out += merge_up_rec(btree, next_state, n_idx);
        out += is_mm_adj(btree, next_state, n_idx) ? loghalf : 0.0; //if we merge just one reverse prob is different as move is always reverse move

    } else if (below_mm && !is_merged)
    {   
        // Node is not merged as a part of an MM & is directly above a multi merger
        // Eligible for single branch merge move
        out += collapse_one(btree, next_state, n_idx);
        out += loghalf; //P choosing reverse move
    }  else if (is_merged && !below_mm) 
    {
        // Node is merged and solely a part of the MM boundary
        // Eligible for branching split move
        out += split_down_rec(btree, next_state, true, n_idx);
    } else if (is_merged) 
    {
        // Node is merged and a part of MM interior
        // Eligible for both single branch split and branching split move
        const double u = std::log(R::runif(0.0,1.0));
        out += -loghalf;

        if (u < loghalf)
        {
            out += split_down_rec(btree, next_state, true, n_idx);
        } else 
        {
            out += expand_one(btree, next_state, n_idx);
        }
    }
    
    times_from_taus(btree, next_state, btree.root_idx());
    
    // Reverse move lp of selecting node
    const double n_lp_rev = tau_lwfunc(next_state, btree, n_idx) - tau_lwsum(next_state, btree);
    out += n_lp_rev;

    return out;
}

double mmctime::mm_prop2::lp_node(const btree& btree, const mc_state& state, int node_idx) const
{
    const double lwsum = tau_lwsum(state, btree);
    const double n_lw = tau_lwfunc(state, btree, node_idx);
    return n_lw-lwsum;
}

double mmctime::mm_prop2::tau_lwfunc(const mc_state& mcs, const btree& btree, int n_idx) const
{
    const int r_idx = btree.root_idx();
    bool can_collapse = (n_idx != r_idx) && 
        (btree.get_parent(n_idx) == btree.root_idx() || _br_data.is_poly(btree.pa_edge(n_idx)));
    // Root should always have prob of zero as no branch above it exists
    return (can_collapse) ? -0.5 * mcs.tau_at(n_idx) * mcs.tau_at(n_idx)  / _sd / _sd : -INFINITY; 
}

double mmctime::mm_prop2::tau_lwsum(const mc_state& mcs, const btree& btree) const
{
    double lwtot = -INFINITY;
    for (int i = btree.n_tip(); i < btree.n_node(); i++)
    {
        lwtot = log_sum_exp(lwtot, tau_lwfunc(mcs, btree, i));
    }
    return lwtot;
}

bool mmctime::mm_prop2::is_mm_adj(const btree& btree,
    const mc_state& new_state,
    int node_idx) const
{
    const int pa = btree.get_parent(node_idx);
    return ((pa != btree.root_idx()) && (new_state.q_at(pa) == 1));
}

std::tuple<double, int> mmctime::mm_prop2::sample_tau_idx(const mc_state& mcs, const btree& btree) const 
{
    double lwtot = tau_lwsum(mcs, btree);

    const double u = std::log(R::runif(0.0,1.0));
    int k = btree.n_tip();
    
    double ll = tau_lwfunc(mcs, btree, k) - lwtot;
    double rsum = ll; 

    while (rsum < u && k < btree.n_node())
    {
        k++; 
        ll = tau_lwfunc(mcs, btree, k) - lwtot;
        rsum = log_sum_exp(rsum, ll);
    }
    return std::make_tuple(ll, k);
}

double mmctime::mm_prop2::merge_up_rec(const btree& btree,
    mc_state& new_state,
    int node_idx) const
{

    double out = 0.0;
    out += collapse_one(btree, new_state, node_idx);
    const int pa = btree.get_parent(node_idx);

    if (!(pa==btree.root_idx()) && !is_mm_adj(btree, new_state, pa))
    {
        const double u = std::log(R::runif(0.0,1.0));
        const double lw = tau_lwfunc(new_state, btree, pa);
        if (u < lw)
        {
            out += -lw;
            out += merge_up_rec(btree, new_state, pa);
        } else
        {
            out += -log_1m_exp(lw);
        }
    }
    return out;
}

double mmctime::mm_prop2::split_down_rec(const btree& btree,
    mc_state& new_state,
    bool first_node,
    int node_idx) const
{
    double out = 0.0;

    out += expand_one(btree, new_state, node_idx);
    if (!first_node)
    {
        const double lw = tau_lwfunc(new_state, btree, node_idx);
        out += lw; // Lp of reverse move including this node
    }

    const int pa = btree.get_parent(node_idx);
    const bool can_merge = !(pa==btree.root_idx()) && !is_mm_adj(btree, new_state, pa);

    if (new_state.q_at(pa)==1)
    {
        out += split_down_rec(btree, new_state, false, pa);
    } else if (can_merge)
    {
        // Else child is expanded
        const double lw = tau_lwfunc(new_state, btree, pa);
        out += log_1m_exp(lw); // Lp of reverse move not including this child
    }

    return out;
}

double mmctime::mm_prop2::expand_one(const btree& btree,
    mc_state& new_state,
    int node_idx) const
{
    double out = 0.0;

    const double tau_prop = sample_tau();
    out += -prop_lp_tau(tau_prop); //Lp of time proposal 
        
    new_state.tau_at(node_idx) = tau_prop;
    new_state.q_at(node_idx) = 0;
    
    return out;
}


double mmctime::mm_prop2::collapse_one(const btree& btree,
    mc_state& new_state,
    int node_idx) const
{
    double out = 0.0;

    const double tau_old = new_state.tau_at(node_idx);
    new_state.q_at(node_idx) = 1;
    new_state.tau_at(node_idx) = 0.0;

    out += prop_lp_tau(tau_old); //Lp of reverse move proposing current value of tau
    
    return out;
}

double mmctime::mm_prop2::sample_tau() const
{
    return std::fabs(R::rnorm(0.0, _sd));
}
double mmctime::mm_prop2::prop_lp_tau(double tau) const
{
    return(std::log(2) + R::dnorm(tau, 0.0, _sd, true));
}
