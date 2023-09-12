#ifndef PROB_UTILS_H // include guard
#define PROB_UTILS_H

#include <cmath>
#include <stdexcept>
#include <algorithm>
namespace mmctime
{
    double log_1m_exp(double x);
    double log_sum_exp(double a, double b);

    template<class T>
    double log_sum_exp(const T& vec); 
}

template<class T>
double mmctime::log_sum_exp(const T& vec)
{
    if(vec.size() == 0)
        std::logic_error("Vector must be non-empty");
    double a = *std::max_element(vec.cbegin(), vec.cend());
    double s = 0.0;
    for (auto it = vec.cbegin(); it != vec.cend(); ++it) 
    {
        s += std::exp(*it - a);
    }
    return a + std::log(s);
}

#endif