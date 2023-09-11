#include "coal_models.h"
#include "probs.h"
#include "prob_utils.h"
#include "state_manip.h"
#include "lambda.h"
#include <tuple>

typedef mmctime::km_beta_model kmb;
typedef mmctime::ds_model ds;
typedef mmctime::beta_model beta;
typedef mmctime::km_model km;

inline double exp_transf(double x) {return std::exp(x);}
inline double logis_transf(double x) {return 1.0 / (1.0 + std::exp(-x));}

inline double jexp(double x) {return x;}
inline double jlogit(double x) {return x - 2.0 * mmctime::log_sum_exp(0.0, x);}

inline double count_mm(const mmctime::btree& btree, const mmctime::mc_state& mc_state)
{
    auto event_rep = convert_state(mc_state, btree);
    int count = 0;
    for (auto e : event_rep)
    {
        count += std::get<1>(e) < -1 ? 1 : 0;
    }
    return count;
}

inline double max_msize(const mmctime::btree& btree, const mmctime::mc_state& mc_state)
{
    auto event_rep = convert_state(mc_state, btree);
    int max_msize = 0;
    for (auto e : event_rep)
    {
        if(max_msize < -std::get<1>(e)) max_msize = -std::get<1>(e);
    }
    return max_msize + 1;
}

/*
KM_BETA MODEL
*/
double kmb::transform_pars(par_state ps) const
{
    ps.constr_at(0) = exp_transf(ps.unconstr_at(0)); //nu
    ps.constr_at(1) = logis_transf(ps.unconstr_at(1)); //phi
    ps.constr_at(2) = logis_transf(ps.unconstr_at(2)); //alpha

    return jexp(ps.unconstr_at(0)) + jlogit(ps.unconstr_at(1)) + jlogit(ps.unconstr_at(2));
}

double kmb::prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    const double nu = ps.constr_at(0); 
    const double phi = ps.constr_at(1); 
    const double alpha = ps.constr_at(2);

    const km_beta_lambda lambda{phi, nu, alpha+1.0};
    return prior_lp_lambda(btree, mc_state, lambda);
}

std::vector<double> kmb::summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const
{   
    return std::vector<double>{count_mm(btree, mc_state), max_msize(btree, mc_state)};
}

/*
DURRET-SCHWEINSBERG MODEL
*/
double ds::transform_pars(par_state ps) const 
{
    ps.constr_at(0) = exp_transf(ps.unconstr_at(0)); //nu
    ps.constr_at(1) = logis_transf(ps.unconstr_at(1)); //phi

    return jexp(ps.unconstr_at(0)) + jlogit(ps.unconstr_at(1));
}

double ds::prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    const double nu = ps.constr_at(0); 
    const double phi = ps.constr_at(1); 

    const ds_lambda lambda{phi, nu};

    return prior_lp_lambda(btree, mc_state, lambda);
}

std::vector<double> ds::summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    return std::vector<double>{count_mm(btree, mc_state), max_msize(btree, mc_state)};
}
/*
BETA MODEL
*/
double beta::transform_pars(par_state ps) const
{
    ps.constr_at(0) = exp_transf(ps.unconstr_at(0)); //nu
    ps.constr_at(1) = logis_transf(ps.unconstr_at(1)); //alpha

    return jexp(ps.unconstr_at(0)) + jlogit(ps.unconstr_at(1));
}

double beta::prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    const double nu = ps.constr_at(0); 
    const double alpha = ps.constr_at(1);

    const beta_lambda lambda{nu, 1.0+alpha};
    return prior_lp_lambda(btree, mc_state, lambda);
}

std::vector<double> beta::summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    return std::vector<double>{count_mm(btree, mc_state), max_msize(btree, mc_state)};
}
/*
KINGMAN MODEL
*/
double km::transform_pars(par_state ps) const
{
    ps.constr_at(0) = exp_transf(ps.unconstr_at(0)); //nu
    return jexp(ps.unconstr_at(0));
}

double km::prior_lp_tree(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    const double nu = ps.constr_at(0); 
    const kingman coal{nu};
    return prior_lp_lambda(btree, mc_state, coal);
}

std::vector<double> km::summaries(const btree& btree, const mc_state& mc_state, const par_state ps) const
{
    return std::vector<double>();
}