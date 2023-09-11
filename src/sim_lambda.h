#ifndef SIM_LAMBDA
#define SIM_LAMBDA

#include "lambda.h"
#include "prob_utils.h"
#include <Rcpp.h>
#include <cmath>
#include <stdexcept>

namespace mmctime
{
    template <typename T>
    Rcpp::List sim_lambda_direct(Rcpp::NumericVector samp_times, Rcpp::IntegerVector n_samp, const T& coal);
}

template <typename T>
Rcpp::List mmctime::sim_lambda_direct(Rcpp::NumericVector samp_times, Rcpp::IntegerVector n_samp, const T& coal)
{
    int i = 0;
    const int m = n_samp.size();
    double t = samp_times(i);
    int at = n_samp(i);

    i++;

    Rcpp::NumericVector coal_times;
    Rcpp::IntegerVector merger_sizes; 

    double lp = 0.0;

    while (i < m || at > 1)
    {
        if (at == 1)
        {
            t = samp_times(i);
            at += n_samp(i);
            i++;
        } else 
        {
            const double lrate_tot = log_total_rate(at, coal);
            const double wt = R::rexp(1.0/std::exp(lrate_tot));
            if (i < m && t + wt > samp_times(i))
            {
                const double t_diff = samp_times(i)-t;
                lp += R::dpois(0, t_diff * std::exp(lrate_tot), true); //lp P[wt+t > next sampling time]

                t = samp_times(i);
                at += n_samp(i);
                i++;
            } else 
            {
                lp += R::dexp(wt, 1.0/std::exp(lrate_tot), true); //lp P[wt = x]
                t += wt;
                const double u = std::log(R::runif(0.0,1.0));
                int k = 1;
                double rsum = -INFINITY;
                double log_binom = std::log(at); 
                double ll = 0.0; 
                double max_ll = -INFINITY;
                while (rsum < u && k < at)
                {
                    k++;
                    log_binom += std::log(at - k + 1) - std::log(k);
                    ll = coal.log_lambda_bk(at,k) + log_binom - lrate_tot;
                    max_ll = std::max(max_ll, ll);
                    rsum = log_sum_exp(rsum, ll);
                }
                if(ll == -INFINITY) 
                    throw std::runtime_error("Critical floating point underflow in simulator! Total rate: " + 
                        std::to_string(rsum) + 
                        "Max event rate: " + 
                        std::to_string(max_ll) + 
                        "u: " +
                        std::to_string(u));
                lp += ll; //lp P[k merger | wt]
                lp -= log_binom; //lp P[which nodes merge | k merger]
                coal_times.push_back(t);
                merger_sizes.push_back(k);
                at -= (k-1);
            }
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("coal_times")=coal_times, 
        Rcpp::Named("merger_sizes")=merger_sizes,
        Rcpp::Named("sim_lp")=lp);
}

#endif