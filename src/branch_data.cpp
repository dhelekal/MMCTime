#include "branch_data.h"

int mmctime::branch_data_IS::n_branch() const
{
    return _branch_vec.size();
}

int mmctime::branch_data_IS::get_branch_muts(int branch) const
{
    return _branch_vec.at(branch);
}

bool mmctime::branch_data_IS::is_poly(int branch, double tol) const
{
    return _branch_vec.at(branch) < tol;
}
