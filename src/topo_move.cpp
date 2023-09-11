#include "topo_move.h"       
#include "prob_utils.h"
#include "state_manip.h" 
#include <cmath>
#include <vector>
#include <stdexcept>
        
double mmctime::topo_prop::prop(const btree& curr_btree, const mc_state& curr_state, btree& next_btree, mc_state& next_state) const
{
    double qr = 0.0;
    const int n = pivot_set.size();

    next_state.copy_from(curr_state);
    next_btree.copy_from(curr_btree);

    const double u1 = R::runif(0.0,1.0);

    const int pivot = pivot_set.at(std::trunc(n * u1));
    if (!br_data.is_poly(next_btree.pa_edge(pivot))) throw 
        std::logic_error("Pivot must be a polytomy");

    const int pa_n = next_btree.get_parent(pivot);

    int pa_c1, pa_c2;
    std::tie(pa_c1, pa_c2) = next_btree.get_children(pa_n);

    const bool is_left = pivot == pa_c1;
    const int other = is_left ? pa_c2 : pa_c1;

    std::vector<int> desc_set_pivot;
    desc_set_pivot.reserve(100);
    _get_poly_desc(pivot, desc_set_pivot, next_btree, false);

    std::vector<int> desc_set_other;
    desc_set_other.reserve(100);
    _get_poly_desc(other, desc_set_other, next_btree, true);

    const int sz_pset = desc_set_pivot.size();
    const int i1 = std::trunc(R::runif(0.0,1.0)*sz_pset);
    const int a = desc_set_pivot.at(i1);

    const int sz_oset = desc_set_other.size();
    const int i2 = std::trunc(R::runif(0.0,1.0)*sz_oset);
    const int b = desc_set_other.at(i2);

    const double q = -std::log(sz_pset)-std::log(sz_oset);
    const double prop_lp = (i2 != 0) ? log_sum_exp(q, -std::log(sz_pset+1)-std::log(sz_oset-1)) : q; //If i2 != 0 could also pick sibling as pivot and obtain same move

    qr += -prop_lp;
    
    next_btree.polytomy_swap(a, b, pa_n);
    times_from_taus(next_btree, next_state, pa_n);

    int pa_c1n, pa_c2n;
    std::tie(pa_c1n, pa_c2n) = next_btree.get_children(pa_n);
    
    const int sz_pnew = _count_poly_desc(is_left ? pa_c1n : pa_c2n, next_btree, false);
    const int sz_onew = _count_poly_desc(is_left ? pa_c2n : pa_c1n, next_btree, true);

    const double r = -std::log(sz_pnew)-std::log(sz_onew);
    const double rev_lp = (i2 != 0) ? log_sum_exp(r, -std::log(sz_pnew+1)-std::log(sz_onew-1)) : r;

    qr += rev_lp;
    return qr;
}

void mmctime::topo_prop::_get_poly_desc(int node, std::vector<int>& v, const btree& btree, bool count_this) const
{    
    if(count_this) v.push_back(node);
    if(!btree.is_tip(node) && br_data.is_poly(btree.pa_edge(node)))
    {
        int n_c1, n_c2;
        std::tie(n_c1, n_c2) = btree.get_children(node);
        
        _get_poly_desc(n_c1, v, btree, true);
        _get_poly_desc(n_c2, v, btree, true);
    }
}

int mmctime::topo_prop::_count_poly_desc(int node, const btree& btree, bool count_this) const
{
    int count = count_this ? 1 : 0;
    if(!btree.is_tip(node) && br_data.is_poly(btree.pa_edge(node)))
    {
        int n_c1, n_c2;
        std::tie(n_c1, n_c2) = btree.get_children(node);
        count += _count_poly_desc(n_c1, btree, true);
        count += _count_poly_desc(n_c2, btree, true);   
    }
    return count;
}
