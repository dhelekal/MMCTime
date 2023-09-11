//lambda.cpp
#include "lambda.h"
#include <cmath>
#include <iostream>

using namespace std;

double mmctime::km_beta_lambda::log_lambda_bk(int b, int k) const
{
    const double lphi = log(phi);
    const double mlphi = log_1m_exp(lphi);
    const double lnu = log(nu);

    const double lbrate = lphi + lgamma(k-alpha) + lgamma(b-k+alpha) - lgamma(b) - lgamma(2.0-alpha) - lgamma(alpha);
    const double lr = ((k==2) ? log_sum_exp(mlphi, lbrate) : lbrate);
    return lr + lnu; 
}

double mmctime::beta_lambda::log_lambda_bk(int b, int k) const
{
    const double lnu = log(nu);
    const double lr = lgamma(k-alpha) + lgamma(b-k+alpha) - lgamma(b) - lgamma(2.0-alpha) - lgamma(alpha);
    return lr + lnu; 
}

double mmctime::bs_lambda::log_lambda_bk(int b, int k) const
{
    const double lnu = log(nu);
    return (lgamma(k-2+1) + lgamma(b-k+1) - lgamma(b)) + lnu;
}

double mmctime::kingman::log_lambda_bk(int b, int k) const
{
    const double lnu = log(nu);
    return k==2 ? lnu : -INFINITY;
}

double mmctime::ds_lambda::log_lambda_bk(int b, int k) const
{
    const double lphi = log(phi);
    const double lnu = log(nu);
    const double mlphi = log_1m_exp(lphi);

    const double lsrate = lgamma(1.0+b-k) + lgamma(k) - lgamma(1.0+b) - log(2.0) + lphi;
    const double lr = ((k==2) ? log_sum_exp(mlphi, lsrate) : lsrate);

    return lr + lnu;
}

double mmctime::orthant_correction(const std::vector<std::tuple<double, int>>& events)
{
    const int n = events.size();
    auto lnbtopo = [](int a) 
    {
        return (std::lgamma(2.0*a-3.0+1.0) - std::lgamma(a-2.0+1.0) - (a-2.0) * std::log(2.0));
    };

    double out = 0.0;
    for (int i = 0; i < n; i++)
    {
        const int n_merge = -std::get<1>(events.at(i)) + 1;
        if (n_merge > 0)
        {
            out += -lnbtopo(n_merge);
        } 
    }
    return out;
}

template<>
double mmctime::lambda_lp(const std::vector<std::tuple<double, int>>& events, int n_tip, kingman lambda)
{
    double lp = 0.0;
    double t_curr = 0.0;

    int at = 0; 

    for (auto e : events)
    {
        int da;
        double t_next;

        std::tie(t_next, da) = e;

        const double r_tot = (at > 1) ? lambda.log_lambda_bk(2,2) + std::log(at) + std::log(at-1) - std::log(2): -INFINITY;

        const double wt = t_next-t_curr;
        const double event_lp = -std::exp(r_tot+std::log(wt)) + (da < 0 ? lambda.log_lambda_bk(at, -da+1) : 0.0);
        lp += event_lp;

        t_curr=t_next;
        at += da;
    }
    if (at != 1) throw std::logic_error("Illegal state reached, A(t) != 1! A(t): " + std::to_string(at));
    return lp;
}