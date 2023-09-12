#ifndef TREE_INIT
#define TREE_INIT

#include <Rcpp.h>

namespace mmctime
{    
    void init_times(const Rcpp::IntegerMatrix& topo_mat,
        const Rcpp::IntegerVector& branch_data,
        int n_node,
        int n_tip,
        double mu,
        int root_br_pos,
        int root_idx,
        Rcpp::NumericVector& node_times);
    
    double init_rec(const Rcpp::IntegerMatrix& topo_mat,
        const Rcpp::IntegerVector& branch_data, 
        const std::vector<std::vector<int> >& ie_vec,
        double mu, 
        int n_tip, 
        int root_br_pos,
        int root_idx,
        Rcpp::NumericVector& node_times, 
        int node_idx);

    Rcpp::IntegerMatrix build_topo_mat(const Rcpp::IntegerMatrix& br_mat, 
        int n_node, 
        int root_pos);

    Rcpp::IntegerMatrix build_branch_mat(const Rcpp::IntegerMatrix& br_mat, 
        int n_node, 
        int root_pos);

    Rcpp::NumericVector create_sbounds(const Rcpp::IntegerMatrix& topo_mat, 
        const Rcpp::NumericVector& tip_times,
        int n_tip,
        int root_idx);
    
    double sbounds_rec(const Rcpp::IntegerMatrix& topo_mat, 
        const Rcpp::NumericVector& tip_times,
        int n_tip,
        Rcpp::NumericVector& sbounds,
        int node_idx);
}

#endif