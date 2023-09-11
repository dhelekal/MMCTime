#include "tree_init.h"
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <functional>

using namespace Rcpp;

void mmctime::init_times(const Rcpp::IntegerMatrix& topo_mat,
    const Rcpp::IntegerVector& branch_data,
    int n_node,
    int n_tip,
    double mu,
    int root_br_pos,
    int root_idx,
    Rcpp::NumericVector& node_times)
{
    std::vector<std::vector<int> > ie_vec(n_node, std::vector<int>());

    for (int i = 0; i < topo_mat.nrow(); i++)
    {
        if (i!=root_idx && (root_br_pos == NA_INTEGER || topo_mat.at(i, 0)-1 != root_idx))
        {
            const int pe = topo_mat.at(i, 3) - 1;
            ie_vec.at(i).push_back(pe);
            ie_vec.at(topo_mat.at(i, 0) - 1).push_back(pe);
        }
    }

    init_rec(topo_mat, 
        branch_data, 
        ie_vec, 
        mu, 
        n_tip, 
        root_br_pos,
        root_idx,
        node_times, 
        root_idx);
}

double mmctime::init_rec(const Rcpp::IntegerMatrix& topo_mat, 
    const Rcpp::IntegerVector& branch_data,
    const std::vector<std::vector<int> >& ie_vec,
    double mu, 
    int n_tip, 
    int root_br_pos,
    int root_idx,
    Rcpp::NumericVector& node_times, 
    int node_idx)
{
    const double eps = 0.01;
    if (node_idx >= n_tip) 
    {   
        const int n_c1 = topo_mat.at(node_idx, 1) - 1;
        const int n_c2 = topo_mat.at(node_idx, 2) - 1;

        const double c1_time = init_rec(topo_mat, branch_data, ie_vec, mu, n_tip, root_br_pos, root_idx, node_times, n_c1);
        const double c2_time = init_rec(topo_mat, branch_data, ie_vec, mu, n_tip, root_br_pos, root_idx, node_times, n_c2);
        
        int mut_1 = NA_INTEGER, mut_2 = NA_INTEGER; 
        int ctr = 0;
        if (root_br_pos == NA_INTEGER || node_idx != root_idx) 
        {
            for (int e : ie_vec.at(node_idx))
            {
                const int pa_e = topo_mat.at(node_idx, 3) - 1;
                if (e != pa_e && ctr == 0)
                {
                    mut_1 = branch_data.at(e);
                    ctr++;
                } else if (e != pa_e)
                {
                    mut_2 = branch_data.at(e);
                }
            }
        } else
        {
            mut_1 = branch_data.at(root_br_pos-1)/2.0;
            mut_2 = branch_data.at(root_br_pos-1)/2.0;
        }

        const double out = std::max(c1_time + std::max(mut_1/mu, eps/mu), c2_time + std::max(mut_2/mu, eps/mu));

        node_times.at(node_idx) = out;
        return out;
    } else
    {
        return node_times.at(node_idx);
    }
}

Rcpp::IntegerMatrix mmctime::build_topo_mat(const Rcpp::IntegerMatrix& br_mat, 
    int n_node, 
    int root_pos,
    double cutoff)
{
    Rcpp::IntegerMatrix topo_mat(n_node, 4);
    std::fill(topo_mat.begin(), topo_mat.end(), NA_INTEGER);

    std::vector<std::vector<int> > ie_vec(n_node, std::vector<int>());
    for (int i = 0; i < br_mat.nrow(); i++)
    {
        ie_vec.at(br_mat.at(i, 0)-1).push_back(i);
        ie_vec.at(br_mat.at(i, 1)-1).push_back(i);
    }

    const int r_pos_ofs = root_pos - 1;
    const int root_idx = n_node - 1;
    
    int rc1;
    int rc2;
    if (root_pos == NA_INTEGER)
    {
        bool first_set = false;
        for (int e : ie_vec.at(root_idx))
        {
            const int e1 = br_mat.at(e, 0) - 1;
            const int e2 = br_mat.at(e, 1) - 1;
            const int other = (e1 != root_idx) ? e1 : e2;
            if (!first_set)
            {
                rc1 = other;
                first_set = true;
            } else
            {
                rc2 = other;
            }
        }
    } else
    {
        rc1 = br_mat.at(r_pos_ofs, 0) - 1;
        rc2 = br_mat.at(r_pos_ofs, 1) - 1;
    }

    topo_mat.at(root_idx, 1) = rc1 + 1;
    topo_mat.at(root_idx, 2) = rc2 + 1;

    std::function<void(int, int)> fill_node = [&ie_vec, &br_mat, &topo_mat, &fill_node, r_pos_ofs, cutoff](int node, int pa_node)  
    {
        topo_mat.at(node, 0) = pa_node + 1;

        int c_idx = 1;
        for (int e : ie_vec.at(node))
        {
            if (e == r_pos_ofs) continue;

            const int e1 = br_mat.at(e, 0) - 1;
            const int e2 = br_mat.at(e, 1) - 1;

            const int other = (e1 != node) ? e1 : e2;
            if (other != pa_node)
            {   
                fill_node(other, node);
                topo_mat.at(node, c_idx) = other + 1;
                c_idx++;
            } else 
            {
                topo_mat.at(node, 3) = e + 1;
            }
        }
    };

    fill_node(rc1, root_idx);
    fill_node(rc2, root_idx);

    return topo_mat;
}

Rcpp::IntegerMatrix mmctime::build_branch_mat(const Rcpp::IntegerMatrix& br_mat, 
    int n_node, 
    int root_pos)
{
    Rcpp::IntegerMatrix out(br_mat.nrow(), 2);
    std::fill(out.begin(), out.end(), NA_INTEGER);

    
    std::vector<std::vector<int> > ie_vec(n_node, std::vector<int>());
    for (int i = 0; i < br_mat.nrow(); i++)
    {
        ie_vec.at(br_mat.at(i, 0)-1).push_back(i);
        ie_vec.at(br_mat.at(i, 1)-1).push_back(i);
    }

    std::function<void(int, int)> fill_branch = [&fill_branch, &br_mat, &ie_vec, &out](int branch, int parent_node)
    {
        const int e1 = br_mat.at(branch, 0) - 1;
        const int e2 = br_mat.at(branch, 1) - 1;

        const int child_node = (e1 == parent_node) ? e2 : e1;
        out.at(branch, 1) = child_node + 1; 

        for (auto e : ie_vec.at(child_node))
        {
            if (e != branch)
            {
                fill_branch(e, child_node);
            }
        }

        out.at(branch, 0) = br_mat.at(branch, 2);
    };

    const int rc1 = br_mat.at(root_pos - 1, 0) - 1;
    const int rc2 = br_mat.at(root_pos - 1, 1) - 1;
    out.at(root_pos-1, 0) = br_mat.at(root_pos-1, 2);

    for (auto e : ie_vec.at(rc1))
    {
        if (e != (root_pos-1))
        {
            fill_branch(e, rc1);
        }
    }

    for (auto e : ie_vec.at(rc2))
    {
        if (e != (root_pos-1))
        {
            fill_branch(e, rc2);
        }
    }
    return out;
}

Rcpp::NumericVector mmctime::create_sbounds(const Rcpp::IntegerMatrix& topo_mat, 
    const Rcpp::NumericVector& tip_times,
    int n_tip,
    int root_idx)
{
    const int n = topo_mat.rows();
    Rcpp::NumericVector out(n, NA_REAL);
    sbounds_rec(topo_mat, tip_times, n_tip, out, root_idx);
    return out;
}

double mmctime::sbounds_rec(const Rcpp::IntegerMatrix& topo_mat, 
    const Rcpp::NumericVector& tip_times,
    int n_tip,
    Rcpp::NumericVector& sbounds,
    int node_idx)
{
    double out;
    if(node_idx < n_tip)
    {
        out = tip_times.at(node_idx);
    } else 
    {
        const int n_c1 = topo_mat.at(node_idx, 1) - 1;
        const int n_c2 = topo_mat.at(node_idx, 2) - 1;
        const double st1 = sbounds_rec(topo_mat, tip_times, n_tip, sbounds, n_c1);
        const double st2 = sbounds_rec(topo_mat, tip_times, n_tip, sbounds, n_c2);

        out = std::max(st1, st2);
    }
    
    sbounds.at(node_idx) = out; 
    return out;
}