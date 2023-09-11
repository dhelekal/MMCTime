#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstring>
#include <iostream>
#include "consense.h"

// [[Rcpp::export]]
Rcpp::IntegerMatrix maj_edge(Rcpp::List trees)
{
    const int n_trees = trees.size();
    const int cutoff = n_trees/2 + 1;
    const mmcutils::clade_table ct(trees);
    const auto mt = mmcutils::make_m_tree(ct, cutoff);

    Rcpp::IntegerMatrix edge = mt.as_edge();
    return edge;
}


// [[Rcpp::export]]
Rcpp::IntegerVector nodes_in_clades(Rcpp::List trees, Rcpp::IntegerMatrix clades)
{
    return mmcutils::count_nodes_in_polys(trees, clades);
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix find_clade_indices(Rcpp::List trees, Rcpp::IntegerMatrix clades)
{
    return mmcutils::find_clade_index(trees, clades);
}

//' @export
// [[Rcpp::export]]
void median_clade_times(Rcpp::List x, Rcpp::List all_trees, Rcpp::NumericMatrix ts)
{
    mmcutils::median_clade_times(x, all_trees, ts);
}