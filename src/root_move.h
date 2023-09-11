
#ifndef ROOT_MOVE
#define ROOT_MOVE

#include "tree_data.h"
#include "par_data.h"

namespace mmctime
{
    class root_prop
    {   
    public:
        double prop(const btree& curr_btree, const mc_state& curr_state, btree& next_btree, mc_state& next_state) const;        
    };
}

#endif