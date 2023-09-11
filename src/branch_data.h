#ifndef BRANCH_DATA
#define BRANCH_DATA

#include <Rcpp.h>

namespace mmctime
{
    class branch_data_IS
    {
    public:
        branch_data_IS(Rcpp::IntegerVector branch_vec) : 
        _branch_vec(branch_vec){};

        int n_branch() const;
        int get_branch_muts(int branch) const;
        bool is_poly(int branch, double tol=1e-8) const;
    private:
        Rcpp::IntegerVector _branch_vec;
    };
}
#endif