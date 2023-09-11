#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <map>

namespace mmcutils
{
    struct rtree
    {
        rtree() = default;
        rtree(Rcpp::List phy);
        Rcpp::IntegerMatrix edge;
        std::vector<std::vector<int>> ie_vec;
        int n_tip=0;
        int n_node=0;
    };
 
    struct clade
    {
        clade() = default;
        clade(std::vector<bool> vec);
        bool is_compatible(const clade& other) const;
        bool contains(const clade& other) const;
        bool is_equal(const clade& other) const;
        int n_tip=0;
        int sz=0;
        std::vector<bool> data;
        std::vector<int> present;
    };

    class clade_table
    {
    public:
        clade_table() = default;
        clade_table(Rcpp::List trees);
        std::vector<clade> clades; // Clades
        std::vector<std::vector<int>> m_sizes; // Merger sizes
        std::vector<int> counts; // Occurrence counts
        std::vector<int> order; //Order of clades by occurrence, descending
        int n_clades=0;
        int n_tip=0;
        int n_draws=0;
    private:
        std::vector<bool> _register_clades(const int idx, const rtree& phy, std::unordered_map<std::vector<bool>, int>& lookup);
    };
    struct clade_node
    {
        int parent = -1;
        int target_size = 2;
        std::vector<int> children;
    };

    class clade_tree
    {
    public:
        clade_tree() = default;
        clade_tree(std::vector<clade> cv, std::vector<int> target_sizes, int n_tip);
        Rcpp::IntegerMatrix as_edge() const;
        bool try_insert(const clade cl, int target_size); 
        
        bool is_resolved() const;
        int unresolved() const;
        void print() const;

        int n_tip = 0;
    private:
        bool _insert_rec(int n_idx, const clade& cl, int target_size, bool check_msize);
        int _unresolved = 0;
        int _root_idx = 0;
        std::vector<clade_node> _nodes;
        std::vector<clade> _clades;
    };
    
    clade_tree make_m_tree(const clade_table& table, int cutoff);

    std::vector<bool> find_clades(int idx, std::vector<std::tuple<int, clade>>& cv, const rtree& tr);
    
    Rcpp::IntegerVector count_nodes_in_polys(Rcpp::List tree, Rcpp::IntegerMatrix clades);
    
    Rcpp::IntegerMatrix find_clade_index(Rcpp::List trees, Rcpp::IntegerMatrix clades);
    void median_clade_times(Rcpp::List x, Rcpp::List all_trees, Rcpp::NumericMatrix ts);
    
    int count_internal(int idx, const rtree& tr);

    //Rcpp::NumericMatrix comp_mm_mat(Rcpp::List trees);

}