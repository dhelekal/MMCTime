
#include <Rcpp.h>
#include <tuple>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <iostream>

#include "tree_data.h"
typedef mmctime::btree bt;

bt bt::from_list(Rcpp::List& list)
{
    const int n_node = list["n_node"];
    const int root_idx = n_node-1;
    return btree(list["topo_mat"], list["sbounds"], list["root_br_pos"], list["n_tip"], root_idx, list["root_fixed"]);
}

std::tuple<int, int> bt::get_children(int node_idx) const
{
    const int n_c1 = _topo_mat.at(node_idx, 1) - 1;
    const int n_c2 = _topo_mat.at(node_idx, 2) - 1;

    return std::make_tuple(n_c1, n_c2);
}

int bt::get_parent(int node_idx) const
{
    return _topo_mat.at(node_idx, 0) - 1;
}

double bt::get_sbound(int node_idx) const
{
    return _sbounds.at(node_idx);
}

int bt::pa_edge(int node_idx) const
{
    return _topo_mat.at(node_idx, 3) - 1;
}

int bt::n_tip() const
{
    return _n_tip;
}

int bt::root_idx() const
{ 
    return _root_idx;
}

int bt::root_br_pos() const
{
     return _root_br_pos.at(0)-1;
}

int bt::n_node() const
{
    return 2 * _n_tip - 1;
}

int bt::root_fixed() const
{
    return _root_fixed;
}

bool bt::is_tip(int node_idx) const
{
    return node_idx < _n_tip;
}

double bt::_upd_sbounds(int node_idx)
{
    if (is_tip(node_idx))
    {
        return _sbounds.at(node_idx);
    } else
    {
        int n_c1, n_c2;
        std::tie(n_c1, n_c2)= get_children(node_idx);

        const double s1 = _upd_sbounds(n_c1);
        const double s2 = _upd_sbounds(n_c2);

        const double sb = std::max(s1, s2);
        _sbounds.at(node_idx) = sb;

        return sb;
    }
}

void bt::_swap_child(int pa, int c_old, int c_new)
{
        const int p_c1 = _topo_mat.at(pa, 1) - 1;
        const int p_c2 = _topo_mat.at(pa, 2) - 1;

        if(p_c1 == c_old)
        {
            _topo_mat.at(pa, 1) = c_new + 1;
        } else if(p_c2 == c_old)
        {
            _topo_mat.at(pa, 2) = c_new + 1;
        } else
        {
            throw std::logic_error("Child not found");
        }
        _topo_mat.at(c_new, 0) = pa + 1;
}


void bt::polytomy_swap(int node_idx1, int node_idx2, int pivot)
{
    const int pa1 = get_parent(node_idx1);
    const int pa2 = get_parent(node_idx2);

    _swap_child(pa1, node_idx1, node_idx2);
    _swap_child(pa2, node_idx2, node_idx1);

    _upd_sbounds(pivot);
}

void bt::root_nni(int node_idx, int pivot)
{
    if(_root_fixed) throw
        std::invalid_argument("Cannot reroot a tree with fixed root!");
    if (is_tip(pivot)) throw
        std::invalid_argument("Pivot cannot be a tip!");
    int r_c1, r_c2;
    std::tie(r_c1, r_c2) = get_children(root_idx());

    if (pivot != r_c1 && pivot != r_c2) throw
        std::invalid_argument("Pivot must be a child of the root!");
    
    int p_c1, p_c2;
    std::tie(p_c1, p_c2) = get_children(pivot);

    if (node_idx != p_c1 && node_idx != p_c2) throw
        std::invalid_argument("Node must be a child of pivot!");
    const int other = (pivot == r_c1) ? r_c2 : r_c1;
    const int root_br_pos_new = pa_edge(node_idx);

    //change neighbour topology
    _swap_child(root_idx(), other, node_idx);
    _swap_child(pivot, node_idx, other);

    //set parent edges
    _topo_mat.at(other, 3) = _root_br_pos.at(0);
    _root_br_pos.at(0) = root_br_pos_new + 1;

    _topo_mat.at(node_idx, 3) = NA_INTEGER;

    int np_c1, np_c2;
    std::tie(np_c1, np_c2) = get_children(pivot);

    const double sb = std::max(_sbounds.at(np_c1), _sbounds.at(np_c2));
    _sbounds.at(pivot) = sb;
}

bt bt::copy() const
{
    return btree(Rcpp::clone(_topo_mat), 
        Rcpp::clone(_sbounds),
        Rcpp::clone(_root_br_pos),
        _n_tip,
        _root_idx,
        _root_fixed);
}

void bt::copy_from(const btree& other)
{
    if(other._topo_mat.size() != _topo_mat.size() || 
        other._sbounds.size() != _sbounds.size() || 
        other._root_br_pos.size() != _root_br_pos.size() || 
        _n_tip != other._n_tip ||
        _root_idx != other._root_idx ||
        _root_fixed != other._root_fixed) throw
        std::invalid_argument("Tree dimensions don't match.");
    
    std::copy(other._root_br_pos.begin(), other._root_br_pos.end(), _root_br_pos.begin());
    std::copy(other._topo_mat.begin(), other._topo_mat.end(), _topo_mat.begin());
    std::copy(other._sbounds.begin(), other._sbounds.end(), _sbounds.begin());
}

Rcpp::List bt::as_list()
{
    return Rcpp::List::create(
        Rcpp::Named("topo_mat") = _topo_mat,
        Rcpp::Named("sbounds") = _sbounds,
        Rcpp::Named("root_br_pos") = _root_br_pos,
        Rcpp::Named("n_tip") = _n_tip,
        Rcpp::Named("n_node") = 2*_n_tip - 1,
        Rcpp::Named("root_fixed") = _root_fixed
    );
}