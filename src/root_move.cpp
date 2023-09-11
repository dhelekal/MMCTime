#include "root_move.h"       
#include "prob_utils.h"
#include "state_manip.h" 

#include <cmath>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <iostream>
#include <Rcpp.h>
        
double mmctime::root_prop::prop(const btree& curr_btree, const mc_state& curr_state, btree& next_btree, mc_state& next_state) const
{
    next_state.copy_from(curr_state);
    next_btree.copy_from(curr_btree);

    const int r_idx = next_btree.root_idx();

    int r_c1, r_c2;
    std::tie(r_c1, r_c2) = next_btree.get_children(r_idx);

    const double lhalf = std::log(0.5);
    double qr = 0.0;
    int pivot;

    if (next_btree.is_tip(r_c1) || next_btree.is_tip(r_c2))
    {
        if (next_btree.is_tip(r_c1) && next_btree.is_tip(r_c2)) throw 
            std::logic_error("Both root descendants are tips!");
        // Only one outgroup can be a tip in trees with more than 2 tips.
        pivot = next_btree.is_tip(r_c1) ? r_c2 : r_c1;
    } else if (next_state.q_at(r_c1) != next_state.q_at(r_c2))
    {
        // If only one branch is collapsed always use it's descendant as pivot
        pivot = (next_state.q_at(r_c1) == 1) ? r_c1 : r_c2;
    } else
    {
        pivot = R::runif(0.0, 1.0) > 0.5 ? r_c1 : r_c2;
        qr += -lhalf;
    }

    // Select which branch to move root to
    int p_c1, p_c2;
    std::tie(p_c1, p_c2) = next_btree.get_children(pivot);
    const int outnode = R::runif(0.0, 1.0) > 0.5 ? p_c1 : p_c2;
    
    //Interchange
    next_btree.root_nni(outnode, pivot);

    //Update prop ratio
    int nr_c1, nr_c2;
    std::tie(nr_c1, nr_c2) = next_btree.get_children(r_idx);

    if((!next_btree.is_tip(nr_c1) && !next_btree.is_tip(nr_c2)) && 
        next_state.q_at(nr_c1) ==  next_state.q_at(nr_c2))
    {
        qr += lhalf;
    }

    times_from_taus(next_btree, next_state, r_idx);
    return qr;
}

