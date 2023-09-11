#include "prob_utils.h"
#include <limits>
#include <cmath>
#include <cstring>

double mmctime::log_1m_exp(double x)
{
    //From `Accurately Computing log(1 − exp(−|a|))` by Martin Machler
    double out = NAN; 
    if (std::abs(x) < 1e-8)
    {
        out = -INFINITY;
    } else
    {
        if (x > 0) throw 
            std::logic_error("X must be negative. X: " + std::to_string(x));
        double a = -std::log(2);
        if (x > a) 
        {
            out = std::log(-std::expm1(x));
        } else 
        {
            out = std::log1p(-std::exp(x)); 
        }
    }
    return out;
}

double mmctime::log_sum_exp(double a, double b)
{
    const double m = std::max(a,b);
    double out; 
    if (m == -INFINITY)
    {
        out = -INFINITY;
    } else
    {
        const double s = std::exp(a-m) + std::exp(b-m);
        out = m + std::log(s);
    }
    return out;
}