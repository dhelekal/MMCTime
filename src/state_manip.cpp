#include "state_manip.h"
#include <cmath>
#include <stdexcept>
#include <vector>

std::vector<std::tuple<double,int>> mmctime::convert_state(const mc_state& mc_state,
    const btree& btree)
{
    const int n = mc_state.n_times();
    
    std::vector<std::tuple<double,int>> out;
    out.reserve(n);

    for (int i = 0; i < n; i++)
    {
        if (mc_state.q_at(i) == 0)
        {
            int da = q_rec(mc_state, btree, i);
            out.push_back(std::make_tuple(mc_state.t_at(i), da));
        }       
    }

    return out;
}
    
int mmctime::q_rec(const mc_state& mc_state,
    const btree& btree,
    int node_idx)
{
    int out;
    if (btree.is_tip(node_idx))
    {
        out = 1;
    } else
    {
        int n_c1, n_c2;
        std::tie(n_c1, n_c2) = btree.get_children(node_idx);

        out = -1;

        if (mc_state.q_at(n_c1) == 1) out += q_rec(mc_state, btree, n_c1);
        if (mc_state.q_at(n_c2) == 1) out += q_rec(mc_state, btree, n_c2);
    }
    return out;
}

void mmctime::mm_update_taus(const btree& btree,
    mc_state& mc_state,
    int node_idx,
    double t)
{
    mc_state.tau_at(node_idx) = 0.0;
    
    int n_c1, n_c2;
    std::tie(n_c1, n_c2) = btree.get_children(node_idx);

    if (mc_state.q_at(n_c1) == 1) mm_update_taus(btree, mc_state, n_c1, t);
    if (mc_state.q_at(n_c2) == 1) mm_update_taus(btree, mc_state, n_c2, t);
}

double mmctime::mm_lower_bound(const mc_state& mc_state,
    const btree& btree,
    int node_idx)
{
    int n_c1, n_c2;
    std::tie(n_c1, n_c2) = btree.get_children(node_idx);

    const double t1 = mc_state.q_at(n_c1) ? mm_lower_bound(mc_state, btree, n_c1) : mc_state.t_at(n_c1);
    const double t2 = mc_state.q_at(n_c2) ? mm_lower_bound(mc_state, btree, n_c2) : mc_state.t_at(n_c2);

    return std::min(t1,t2);
}

int mmctime::mm_count(const mc_state& mc_state,
    const btree& btree,
    int node_idx)
{
    int out = 1;
    int n_c1, n_c2;
    std::tie(n_c1, n_c2) = btree.get_children(node_idx);

    out += mc_state.q_at(n_c1) ? mm_count(mc_state, btree, n_c1) : 0;
    out += mc_state.q_at(n_c2) ? mm_count(mc_state, btree, n_c2) : 0;
    
    return out;
}

void mmctime::taus_from_times(const btree& btree, mc_state& mc_state, int node_idx)
{
    if (node_idx == btree.root_idx())
    {
        const double t = mc_state.t_at(node_idx);
        const double sb = btree.get_sbound(node_idx);
        mc_state.tau_at(node_idx) = std::log(t - sb);
    } else 
    {
        const double pa_t = mc_state.t_at(btree.get_parent(node_idx));
        const double sb = btree.get_sbound(node_idx);
        const double t = mc_state.t_at(node_idx);

        mc_state.tau_at(node_idx) = -std::log((t-sb)/(pa_t-sb));
        //-(std::log(t - (pa_t - sb)*sb) - std::log(pa_t - sb));
    }

    int n_c1, n_c2;
    std::tie(n_c1, n_c2) = btree.get_children(node_idx);
    
    if (!btree.is_tip(n_c1)) taus_from_times(btree, mc_state, n_c1);
    if (!btree.is_tip(n_c2)) taus_from_times(btree, mc_state, n_c2);
}

void mmctime::times_from_taus(const btree& btree, mc_state& mc_state, int node_idx)
{
    if (node_idx == btree.root_idx())
    {
        const double tau = mc_state.tau_at(node_idx);
        const double sb = btree.get_sbound(node_idx);
        mc_state.t_at(node_idx) = std::exp(tau) + sb;
    } else 
    {
        const double pa_t = mc_state.t_at(btree.get_parent(node_idx));
        const double sb = btree.get_sbound(node_idx);
        const double tau = mc_state.tau_at(node_idx);

        mc_state.t_at(node_idx) = (pa_t - sb)*std::exp(-tau) + sb;
    }

    int n_c1, n_c2;
    std::tie(n_c1, n_c2) = btree.get_children(node_idx);
    
    if (!btree.is_tip(n_c1)) times_from_taus(btree, mc_state, n_c1);
    if (!btree.is_tip(n_c2)) times_from_taus(btree, mc_state, n_c2);
}

int mmctime::q_count(const mc_state& mc_state)
{
    int out = 0;
    for (int i = 0; i < mc_state.n_times(); i++)
    {
        out += mc_state.q_at(i);
    }
    return out;
}


void mmctime::squash_brs(const btree& btree, double tol, mc_state& mc_state, int node_idx)
{
    if (node_idx != btree.root_idx())
    {
        double t = mc_state.t_at(node_idx);
        double t_pa = mc_state.t_at(btree.get_parent(node_idx));

        if (std::fabs(t_pa-t) <= tol)
        {
                mc_state.q_at(node_idx) = 1;
            mc_state.tau_at(node_idx) = 0.0; 
        }
    }
    int n_c1, n_c2;
    std::tie(n_c1, n_c2) = btree.get_children(node_idx);
        
    if (!btree.is_tip(n_c1)) squash_brs(btree, tol, mc_state, n_c1);
    if (!btree.is_tip(n_c2)) squash_brs(btree, tol, mc_state, n_c2);
}

bool mmctime::validate_times(const btree& btree, const mc_state& mc_state, int node_idx)
{
    const double tol = 1e-8;
    if (btree.is_tip(node_idx))
    {
        return true;
    } else
    {
        int n_c1, n_c2;
        std::tie(n_c1, n_c2) = btree.get_children(node_idx);
        const double t = mc_state.t_at(node_idx);

        bool res1 = t >= mc_state.t_at(n_c1) || std::fabs(t-mc_state.t_at(n_c1)) < tol;
        bool res2 = t >= mc_state.t_at(n_c2) || std::fabs(t-mc_state.t_at(n_c2)) < tol;

        return res1 && res2 && validate_times(btree, mc_state, n_c1) && validate_times(btree, mc_state, n_c2);
    }
}



