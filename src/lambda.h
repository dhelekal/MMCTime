//lambda.h
#ifndef LAMBDA_H // include guard
#define LAMBDA_H

#include "prob_utils.h"
#include <cmath>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <cstring>
#include <iostream>
namespace mmctime
{   
    struct km_beta_lambda
    {
        double phi;
        double nu;
        double alpha;
        double log_lambda_bk(int b, int k) const;
    };

    struct beta_lambda
    {
        double nu;
        double alpha;
        double log_lambda_bk(int b, int k) const;
    };

    struct bs_lambda
    {
        double nu;
        double log_lambda_bk(int b, int k) const;
    };

    struct ds_lambda
    {
        double phi;
        double nu;
        double log_lambda_bk(int b, int k) const; 
    };

    struct kingman
    {
        double nu;
        double log_lambda_bk(int b, int k) const; 
    };

    template <typename T>
    double log_total_rate(int b, const T& coal);

    template<typename T>
    double lambda_lp(const std::vector<std::tuple<double, int>>& events, int n_tip, T lambda);

    double orthant_correction(const std::vector<std::tuple<double, int>>& events);
}

template <typename T>
double mmctime::log_total_rate(int b, const T& coal)
{
    double out = -INFINITY;
    double log_binom = std::log(b); 
    for (int i = 2; i <= b; i++)
    {        
        log_binom = log_binom + std::log(b - i + 1) - std::log(i);
        double lmrate = log_binom + coal.log_lambda_bk(b,i);
        out = log_sum_exp(out, lmrate);
    }
    return out;
}

template<typename T>
double mmctime::lambda_lp(const std::vector<std::tuple<double, int>>& events, int n_tip, T lambda)
{
    double lp = 0.0;
    double t_curr = 0.0;

    int at = 0; 
    
    std::vector<double> rate_vec(n_tip+1, -INFINITY);
    for (int i = 2; i <= n_tip; i++)
    {
        rate_vec.at(i) = log_sum_exp(rate_vec.at(i-1), std::log(i-1) + lambda.log_lambda_bk(i,2));
    }

    for (auto e : events)
    {
        int da;
        double t_next;

        std::tie(t_next, da) = e;

        const double wt = t_next-t_curr;
        const double event_lp = -std::exp(rate_vec.at(at)+std::log(wt)) + (da < 0 ? lambda.log_lambda_bk(at, -da+1) : 0.0);
        lp += event_lp;

        t_curr=t_next;
        at += da;
    }
    if (at != 1) throw std::logic_error("Illegal state reached, at != 1! at: " + std::to_string(at));
    return lp;
}

#endif