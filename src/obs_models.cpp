#include "obs_models.h"
#include <Rcpp.h>
#include "probs.h"

inline double exp_transf(double x) {return std::exp(x);}
inline double jexp(double x) {return x;}

typedef mmctime::arc_model arc;

double arc::transform_pars(par_state ps)
{
    ps.constr_at(0) = exp_transf(ps.unconstr_at(0));
    ps.constr_at(1) = exp_transf(ps.unconstr_at(1));

    return jexp(ps.unconstr_at(0)) + jexp(ps.unconstr_at(1));
}

double arc::obs_lp_branch(const btree& btree, const mc_state& mc_state, const branch_data_IS& obs, const par_state ps) const
{
    struct arc
    {
        double lp(int mut, double len) const
        {
            const double theta = mu / omega; 
            const double p = 1.0/(1.0 + omega);
            const double k = theta * len;
            return R::dnbinom(mut, k, p, true);
        }
        double mu;
        double omega;
    };

    const double mu = ps.constr_at(0);
    const double omega = ps.constr_at(1); 

    arc arc{mu, omega};

    return obs_lp(btree, mc_state, obs, arc);
}