#ifndef TREE_DATA
#define TREE_DATA

#include <Rcpp.h>
#include <tuple>

namespace mmctime
{
    class btree
    {
    public:
        static btree from_list(Rcpp::List& list);
        btree(Rcpp::IntegerMatrix topo_mat, 
            Rcpp::NumericVector sbounds,
            Rcpp::IntegerVector root_br_pos,
            int n_tip, 
            int root_idx,
            int root_fixed) : 
        _topo_mat(topo_mat),  _sbounds(sbounds), _root_br_pos(root_br_pos), _n_tip(n_tip), _root_idx(root_idx), _root_fixed(root_fixed){};
        
        std::tuple<int, int> get_children(int node_idx) const;
        
        int get_parent(int node_idx) const;
        double get_sbound(int node_idx) const;
        int pa_edge(int node_idx) const;

        int root_idx() const;
        int root_br_pos() const;
        bool is_tip(int node_idx) const;
        
        int n_node() const;
        int n_tip() const;
        int root_fixed() const;

        void polytomy_swap(int node_idx1, int node_idx2, int pivot);
        void root_nni(int node_idx, int pivot);
        btree copy() const;
        void copy_from(const btree& other);

        Rcpp::List as_list();
    private:
        double _upd_sbounds(int node_idx);
        void _swap_child(int pa, int c_old, int c_new);

        Rcpp::IntegerMatrix _topo_mat;
        Rcpp::NumericVector _sbounds;
        Rcpp::IntegerVector _root_br_pos;
        int _n_tip;
        int _root_idx;
        int _root_fixed;
    };
}

#endif