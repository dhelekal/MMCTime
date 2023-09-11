
#ifndef TOPO_MOVE
#define TOPO_MOVE

#include "tree_data.h"
#include "par_data.h"
#include "branch_data.h"

#include <Rcpp.h>

namespace mmctime
{
    class topo_prop
    {   
    public:
        Rcpp::IntegerVector pivot_set;
        branch_data_IS br_data;
        double prop(const btree& curr_btree, const mc_state& curr_state, btree& next_btree, mc_state& next_state) const;
    private:
        void _get_poly_desc(int node, std::vector<int>& v, const btree& btree, bool count_this) const;
        int _count_poly_desc(int node, const btree& btree, bool count_this) const;
    };
}

#endif