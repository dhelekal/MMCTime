#include <Rcpp.h>
#include <stdexcept>
#include <cstring>

#include "par_data.h"

typedef mmctime::mc_state mcs;
typedef mmctime::par_state ps;

mcs mcs::from_list(Rcpp::List& list)
{
    return mmctime::mc_state(list["taus"], list["times"], list["qs"], list["n_tip"]);
}

mcs::mc_state(Rcpp::NumericVector taus, Rcpp::NumericVector times, Rcpp::IntegerVector qs, int n_tip) :
    _taus(taus), _times(times), _qs(qs), _sz_times(times.size()), _sz_taus(taus.size()), _n_tip(n_tip)
{

}

int& mcs::q_at(int node_idx)
{
    return _qs.at(node_idx);    
}

double& mcs::tau_at(int node_idx)
{
    return _taus.at(node_idx - _n_tip);
}

double& mcs::t_at(int node_idx)
{
    return _times.at(node_idx);
}

const int& mcs::q_at(int node_idx) const
{
    return _qs.at(node_idx);        
}

const double& mcs::tau_at(int node_idx) const
{
    return _taus.at(node_idx - _n_tip);    
}

const double& mcs::t_at(int node_idx) const
{
    return _times.at(node_idx);
}

int mcs::n_times() const
{
    return _sz_times;
}

int mcs::n_taus() const
{
    return _sz_taus;
}

bool mcs::has_tau(int node_idx) const
{
    return node_idx < _n_tip;
}

mcs mcs::copy() const
{
    return mcs(Rcpp::clone<Rcpp::NumericVector>(_taus), 
        Rcpp::clone<Rcpp::NumericVector>(_times), 
        Rcpp::clone<Rcpp::IntegerVector>(_qs),
        _n_tip);
}

void mcs::copy_from(const mcs& other)
{
    if(other._taus.size() != _taus.size()  || other._times.size() != _times.size() || other._qs.size() != _qs.size() || other._n_tip != _n_tip) throw
        std::invalid_argument("State dimensions don't match.");
    
    std::copy(other._taus.begin(), other._taus.end(), _taus.begin());
    std::copy(other._times.begin(), other._times.end(), _times.begin());
    std::copy(other._qs.begin(), other._qs.end(), _qs.begin());
}

Rcpp::List mcs::as_list()
{
    return Rcpp::List::create(
        Rcpp::Named("taus")=_taus, 
        Rcpp::Named("times")=_times,
        Rcpp::Named("qs")=_qs,
        Rcpp::Named("n_tip")=_n_tip
    );
}

ps ps::from_list(Rcpp::List& list)
{
    return mmctime::par_state(list["pars"], list["transf_pars"]);
}

ps::par_state(Rcpp::NumericVector pars, Rcpp::NumericVector transf_pars) : _pars(pars), _transf_pars(transf_pars), _n_pars(pars.size()), _extent(pars.size())
{
    if(!(_offset + _extent <= _n_pars)) throw
        std::invalid_argument("Invalid view size"); 
}

ps::par_state(Rcpp::NumericVector pars, Rcpp::NumericVector transf_pars, int offset, int extent) : _pars(pars), _transf_pars(transf_pars), _n_pars(pars.size()), _offset(offset), _extent(extent)
{
    if(!(_offset + _extent <= _n_pars)) throw
        std::invalid_argument("Invalid view size"); 
}


int ps::size() const
{
    return _extent;
}

const Rcpp::NumericVector& ps::constr() const
{
    return _transf_pars;
}
Rcpp::NumericVector& ps::constr()
{
    return _transf_pars;
}

const Rcpp::NumericVector& ps::unconstr() const
{
    return _pars;
}
Rcpp::NumericVector& ps::unconstr()
{
    return _pars;
}

const double& ps::constr_at(int i) const
{
    const int idx = i + _offset;
    if (!(idx < _extent+_offset)) throw std::invalid_argument("Index out of bounds! Index: " + 
        std::to_string(i) + " Extent: " + std::to_string(_extent));
    return _transf_pars.at(idx);
}

double& ps::constr_at(int i)
{
    const int idx = i + _offset;
    if (!(idx < _extent+_offset)) throw std::invalid_argument("Index out of bounds! Index: " + 
        std::to_string(i) + " Extent: " + std::to_string(_extent));
    return _transf_pars.at(idx);
}

const double& ps::unconstr_at(int i) const
{
    const int idx = i + _offset;
    if (!(idx < _extent+_offset)) throw std::invalid_argument("Index out of bounds! Index: " + 
        std::to_string(i) + " Extent: " + std::to_string(_extent));
    return _pars.at(idx);
}

double& ps::unconstr_at(int i)
{
    const int idx = i + _offset;
    if (!(idx < _extent+_offset)) throw std::invalid_argument("Index out of bounds! Index: " + 
        std::to_string(i) + " Extent: " + std::to_string(_extent));
    return _pars.at(idx);
}

ps ps::create_view(int offset, int extent)
{
    return par_state(_pars, _transf_pars, offset, extent);
}

const ps ps::create_view(int offset, int extent) const
{
    return par_state(_pars, _transf_pars, offset, extent);
}

ps ps::copy() const
{
    return ps(Rcpp::clone<Rcpp::NumericVector>(_pars), 
        Rcpp::clone<Rcpp::NumericVector>(_transf_pars));
}

void ps::copy_from(const ps& other)
{
    if(other._pars.size() != _pars.size()  || other._transf_pars.size() != _transf_pars.size()) throw
        std::invalid_argument("State dimensions don't match.");
    
    std::copy(other._pars.begin(), other._pars.end(), _pars.begin());
    std::copy(other._transf_pars.begin(), other._transf_pars.end(), _transf_pars.begin());
}

Rcpp::List ps::as_list()
{
    return Rcpp::List::create(
        Rcpp::Named("pars")=_pars, 
        Rcpp::Named("transf_pars")=_transf_pars
    );
}