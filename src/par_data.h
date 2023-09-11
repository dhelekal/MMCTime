#ifndef PAR_DATA
#define PAR_DATA

#include <Rcpp.h>
#include <tuple>

namespace mmctime
{
    class mc_state
    {
    public:
        static mc_state from_list(Rcpp::List& list);
        mc_state(Rcpp::NumericVector taus,
            Rcpp::NumericVector times,
            Rcpp::IntegerVector qs, 
            int n_tip);
        
        int& q_at(int node_idx);
        double& tau_at(int node_idx);
        double& t_at(int node_idx);

        const int& q_at(int node_idx) const;
        const double& tau_at(int node_idx) const;
        const double& t_at(int node_idx) const;

        int n_times() const;
        int n_taus() const;
        bool has_tau(int node_idx) const;

        mc_state copy() const;
        void copy_from(const mc_state& other);

        Rcpp::List as_list();
    private:
        Rcpp::NumericVector _taus;
        Rcpp::NumericVector _times;
        Rcpp::IntegerVector _qs;

        int _sz_times;
        int _sz_taus;
        int _n_tip;
    };
    class par_state
    {
    public:
        static par_state from_list(Rcpp::List& list);
        par_state(Rcpp::NumericVector pars, Rcpp::NumericVector transf_pars);     

        const Rcpp::NumericVector& constr() const;
        Rcpp::NumericVector& constr();

        const Rcpp::NumericVector& unconstr() const;
        Rcpp::NumericVector& unconstr();

        const double& constr_at(int i) const;
        double& constr_at(int i);

        const double& unconstr_at(int i) const;
        double& unconstr_at(int i);

        par_state create_view(int offset, int extent);
        const par_state create_view(int offset, int extent) const;

        int size() const;

        par_state copy() const;
        void copy_from(const par_state& other);

        Rcpp::List as_list();
    private:
        par_state(Rcpp::NumericVector pars, Rcpp::NumericVector transf_pars, int offset, int extent);
        Rcpp::NumericVector _pars;
        Rcpp::NumericVector _transf_pars;
        int _n_pars;
        int _offset = 0;
        int _extent;
    };
}

#endif