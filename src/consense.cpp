#include "consense.h"
#include <tuple>
#include <algorithm>
#include <stdexcept>


mmcutils::clade::clade(std::vector<bool> vec) : n_tip(vec.size()), data(vec)
{
    for (int i = 0; i < data.size(); i++)
    {
        if (vec.at(i)==true) present.push_back(i);
    }    
    sz = present.size();
}

bool mmcutils::clade::is_compatible(const clade& other) const
{
    auto compat = [](const clade& bigger, const clade& smaller)
    {
        bool contains = true;
        bool disjoint = true;
        for (int i = 0; i < smaller.sz; i++)
        {
            contains = contains && bigger.data.at(smaller.present.at(i));
            disjoint = disjoint && !bigger.data.at(smaller.present.at(i));
        }
        return contains || disjoint;
    };

    bool res;
    if (other.sz >= sz)
    {
        res = compat(other, *this);
    } else
    {
        res = compat(*this, other);
    }

    return res;
}

bool mmcutils::clade::contains(const clade& other) const
{
    if (other.sz > sz) return false;
    for (int i = 0; i < other.sz; i++)
    {
        if(!data.at(other.present.at(i))) return false;
    }
    return true;
}

bool mmcutils::clade::is_equal(const clade& other) const
{
    return data == other.data;
}

mmcutils::clade_table::clade_table(Rcpp::List trees)
{
    std::unordered_map<std::vector<bool>, int> lookup;
    n_draws = trees.size();
    for (int i = 0; i < n_draws; i++)
    {
        Rcpp::List tre_list = trees.at(i);
        auto x = rtree(tre_list);

        if (i == 0) n_tip = x.n_tip;
        _register_clades(x.n_tip, x, lookup);
    }

    std::vector<std::tuple<int, int>> occ_idx;
    occ_idx.reserve(n_clades);

    for (int i = 0; i < n_clades; i++)
    {
        occ_idx.push_back(std::make_tuple(i, counts.at(i)));
    }
    std::sort(occ_idx.begin(), occ_idx.end(), [](auto a, auto b) {return std::get<1>(a) > std::get<1>(b);});

    order.reserve(n_clades);
    for (auto x : occ_idx)
    {
        order.push_back(std::get<0>(x));
    }
}

mmcutils::rtree::rtree(Rcpp::List phy)
{ 
    Rcpp::CharacterVector tmp = phy["tip.label"];
    Rcpp::IntegerMatrix edge_tmp = phy["edge"];
    edge = edge_tmp;
    n_tip = tmp.size();
    n_node = edge.nrow()+1;

    ie_vec = std::vector<std::vector<int> >(n_node, std::vector<int>());
    for (int i = 0; i < edge.nrow(); i++)
    {
        ie_vec.at(edge.at(i, 0)-1).push_back(i);
        ie_vec.at(edge.at(i, 1)-1).push_back(i);
    }
}

std::vector<bool> mmcutils::clade_table::_register_clades(const int idx,
    const rtree& phy, 
    std::unordered_map<std::vector<bool>, int>& lookup)
{
    const int n_tip = phy.n_tip;
    std::vector<bool> cl(n_tip, false);
    int msize;
    if (idx < n_tip)
    {
        msize = 1;
        cl.at(idx) = true;
    } else
    {
        msize = 0;
        for (int e : phy.ie_vec.at(idx))
        {
            const int child = phy.edge.at(e, 1) - 1;
            if (child != idx)
            {   
                msize++;
                auto tmp = _register_clades(child, phy, lookup);
                for (int j = 0; j < n_tip; j++)
                {
                    cl.at(j) = cl.at(j) || tmp.at(j);
                }
            } 
        }
    }

    if (lookup.count(cl) == 0)
    {
        const int j = clades.size();
        counts.push_back(1);
        m_sizes.push_back(std::vector<int> {msize});
        clades.push_back(clade(cl));
        lookup.insert({cl, j});
        n_clades++;
    } else
    {
        const int j = lookup.at(cl);
        m_sizes.at(j).push_back(msize);
        counts.at(j)++;
    }
    return cl;
}

mmcutils::clade_tree::clade_tree(std::vector<clade> cv, std::vector<int> target_sizes, int n_tip) : n_tip(n_tip)
{
    const int n = cv.size();
    std::vector<int> ord(n);

    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&cv](auto a, auto b) {return cv.at(a).sz > cv.at(b).sz;});

    _clades.reserve(n); 
    _nodes.reserve(n);

    _root_idx = 0;
    _clades.push_back(cv.at(ord.at(0))); 

    clade_node root_node;
    root_node.target_size = target_sizes.at(ord.at(0));
    _nodes.push_back(root_node);

    for (int i = 1; i < n; i++)
    {
        if(!_insert_rec(_root_idx, cv.at(ord.at(i)), target_sizes.at(ord.at(i)), false)) throw
            std::logic_error("Majority tree not consistent");
    }

    _unresolved = 0;
    for (int i = 0; i < _nodes.size(); i++)
    {
        if (_nodes.at(i).target_size < _nodes.at(i).children.size()) _unresolved++;
    }      
}

bool mmcutils::clade_tree::try_insert(const clade cl, int target_size)
{
    return _insert_rec(_root_idx, cl, target_size, true);
}

bool mmcutils::clade_tree::_insert_rec(int n_idx, const clade& cl, int target_size, bool check_msize)
{
    bool compat = true;
    auto& curr_node = _nodes.at(n_idx);
    const auto& curr_clade = _clades.at(n_idx);
    if (!curr_clade.contains(cl)) throw
        std::logic_error("Illegal state encountered, current clade does not contain candidate.");
    
    for (auto c_idx : curr_node.children)
    {
        const auto& cn = _clades.at(c_idx);
        if (cn.contains(cl))
        {
            //Check if a child contains candidate, attempt to insert it there
            return _insert_rec(c_idx, cl, target_size, check_msize);
        } else
        {
            //Check if candidate clade contradicts any of the existing child nodes
            compat = compat && cn.is_compatible(cl);
        }
    }

    
    const int curr_sz = curr_node.children.size();
    if (compat)
    {
        std::vector<int> candidate_contains;
        std::vector<int> prop_children;
 
        for (int i = 0; i < curr_sz; i++)
        {
            
            const int j = curr_node.children.at(i);
            if(cl.contains(_clades.at(j)))
            {
                candidate_contains.push_back(j);
            } else
            {
                prop_children.push_back(j);
            }
            
        }

        const int prop_sz = prop_children.size() + 1;
        if (!check_msize || prop_sz >= curr_node.target_size)
        {
            
            clade_node tmp; 
            tmp.target_size = target_size;
            tmp.children = candidate_contains;
            tmp.parent = n_idx;
                   
            _clades.push_back(cl);    
            prop_children.push_back(_nodes.size());
            curr_node.children = prop_children;

            if (prop_sz != curr_node.children.size())
                throw std::logic_error("Size doesnt match");

            if (prop_sz == curr_node.target_size)
            {
                _unresolved--;
            }

            if (candidate_contains.size() > target_size)
            {
                _unresolved++;
            }

            _nodes.push_back(tmp);
            return true;
        } else
        {
            std::cout <<  "Tree would be overresolved" << std::endl;
            return false;
        }
    } else
    {
        //Candidate clade contradicts existing tree
        return false;
    }
}\

Rcpp::IntegerMatrix mmcutils::clade_tree::as_edge() const
{
    auto cl_copy = _clades;
    const int n_node = cl_copy.size();
    const int m = n_node-1;

    Rcpp::IntegerMatrix out(m, 2);
    std::fill(out.begin(), out.end(), NA_INTEGER); 
   
    std::sort(cl_copy.begin(), cl_copy.end(), [] (const auto& a, const auto& b) {return a.sz < b.sz;});
    std::sort(cl_copy.begin(), cl_copy.begin() + n_tip, [] (const auto& a, const auto& b) {return a.present.at(0) < b.present.at(0);});
    
    const int r_idx = n_tip;

    for (int i = 0; i < m; i++) //skip root
    {
        const auto& child = cl_copy.at(i);
        for(int j = i + 1; j < n_node; j++)
        {
            if (cl_copy.at(j).contains(child))
            {
                out.at(i, 0) = j + 1;
                out.at(i, 1) = i + 1; 
                break; 
            }
        }
    }
    // Reindex root for ape
    for (int i = 0; i < m; i++) //skip root
    {
        if (out.at(i, 1) == (r_idx + 1))
        {
            out.at(i, 1) = m + 1;
        }

        const int a = out.at(i, 0);
        if (a == (r_idx + 1))
        {
            out.at(i, 0) = m + 1;
        } else if (a == (m + 1))
        {
            out.at(i, 0) = r_idx + 1;
        }
    }
    return out;
}

void mmcutils::clade_tree::print() const
{
    for (auto cl : _clades)
    {
        for (auto p : cl.present)
        {
            std::cout << std::to_string(p); 
            std::cout << " ";
        }
        std::cout << std::endl;
    }
}

bool mmcutils::clade_tree::is_resolved() const
{
    return !(_unresolved > 0);
}

int mmcutils::clade_tree::unresolved() const
{
    return _unresolved;
}

mmcutils::clade_tree mmcutils::make_m_tree(const clade_table& table, int cutoff)
{
    std::vector<clade> tree_clades;
    std::vector<int> target_sizes;
    
    std::vector<int> mm_tmp;
    mm_tmp.reserve(table.n_draws);

    auto mm_med = [&mm_tmp] (const std::vector<int>& v)
    {
        
        mm_tmp = v;
        int k = mm_tmp.size()/2;
        k += (mm_tmp.size() % 2) - 1;
        if (!(k<mm_tmp.size())) throw 
            std::logic_error("Median index exceeds vector size!");
        std::nth_element(mm_tmp.begin(), mm_tmp.begin() + k, mm_tmp.end());
        
        return mm_tmp.at(k);
    };
    for (int i = 0; i < table.n_clades; i++)
    {
        const int j = table.order.at(i);
        if (table.counts.at(j) >= cutoff)
        {
            target_sizes.push_back(mm_med(table.m_sizes.at(j)));
            tree_clades.push_back(table.clades.at(j));
        } else
        {
            break;
        }
    }

    clade_tree tr(tree_clades, target_sizes, table.n_tip);

    return tr;
}


std::vector<bool> mmcutils::find_clades(int idx, std::vector<std::tuple<int, clade>>& cv, const rtree& tr)
{
    const int n_tip = tr.n_tip;
    std::vector<bool> cl(n_tip, 0);
    if (idx < n_tip)
    {
        cl.at(idx) = 1;
    } else
    {
        for (int e : tr.ie_vec.at(idx))
        {
            const int child = tr.edge.at(e, 1) - 1;
            if (child != idx)
            {   
                auto tmp = find_clades(child, cv, tr);
                for (int j = 0; j < n_tip; j++)
                {
                    cl.at(j) = cl.at(j) || tmp.at(j);
                }
            } 
        }
    }
    cv.push_back(std::make_tuple(idx, clade(cl)));
    return cl;
}

Rcpp::IntegerVector mmcutils::count_nodes_in_polys(Rcpp::List tree, Rcpp::IntegerMatrix clades)
{
    
    const int n_clade = clades.nrow();

    std::vector<clade> to_test;
    to_test.reserve(n_clade);

    const int n_tip = clades.ncol();

    for (int i = 0; i < n_clade; i++)
    {
        std::vector<bool> tmp;
        tmp.reserve(n_tip);
        for (int j = 0; j < clades.ncol(); j++)
        {
            const bool val = clades.at(i,j) > 0 ? 1 : 0;
            tmp.push_back(val);
        }
        to_test.push_back(clade(tmp));
    }
    rtree tr(tree);
    const int r_idx = tr.n_tip;
    std::vector<std::tuple<int, clade>> cv;
    cv.reserve(tr.n_node);
    find_clades(r_idx, cv, tr);

    //match clades with nodes and count internal nodes in clade subtrees.
    std::vector<int> node_counts;

    for (const auto& x : to_test)
    {
        for (const auto& y : cv)
        {
            if (x.is_equal(std::get<1>(y)))
            {
                node_counts.push_back(count_internal(std::get<0>(y),tr));
                break;
            }
        }        
    }
    
    std::vector<int> postorder(n_clade);
    std::iota(postorder.begin(), postorder.end(), 0);
    std::sort(postorder.begin(), postorder.end(), [&to_test](auto a, auto b) {return to_test.at(a).sz < to_test.at(b).sz;});

    for (int i = 0; i < n_clade; i++)
    {
        const int ii = postorder.at(i);
        for (int j = 0; j < i; j++)
        {
            const int jj = postorder.at(j);

            if (to_test.at(ii).contains(to_test.at(jj)))
            {
                node_counts.at(ii) -= node_counts.at(jj);
            }
        }
    }
    return Rcpp::wrap(node_counts);
}

int mmcutils::count_internal(int idx, const rtree& tr)
{
    const int n_tip = tr.n_tip;
    int out = 0;
    //if (idx >= n_tip)
    //{
        for (int e : tr.ie_vec.at(idx))
        {
            const int child = tr.edge.at(e, 1) - 1;
            if (child != idx && child >= n_tip)
            {           
                out++;
                out += count_internal(child, tr);
            } 
        }
    //}
    return out;
}

Rcpp::IntegerMatrix mmcutils::find_clade_index(Rcpp::List trees, Rcpp::IntegerMatrix clades)
{
    const int n_clades = clades.nrow();
    const int n_tip = clades.ncol();
    std::vector<clade> to_test;
    to_test.reserve(n_clades);

    std::unordered_map<std::vector<bool>, int> lookup;
 
    for (int i = 0; i < n_clades; i++)
    {
        std::vector<bool> tmp;
        tmp.reserve(n_tip);
        for (int j = 0; j < n_tip; j++)
        {
            const bool val = clades.at(i,j) > 0 ? 1 : 0;
            tmp.push_back(val);
        }
        lookup.insert({tmp, i});
    }

    const int n_trees = trees.length();

    Rcpp::IntegerMatrix out(n_trees, n_clades);
    std::fill(out.begin(), out.end(), NA_INTEGER);
    std::vector<std::tuple<int, clade>> cv;    
    cv.reserve(2*n_tip - 1);

    for (int i = 0; i < n_trees; i++)
    {
        Rcpp::List tre_list = trees.at(i);
        rtree tr(tre_list);
        find_clades(n_tip, cv, tr);
        for (int j = 0; j < cv.size(); j++)
        {
            const auto& tmp = cv.at(j);
            const int a = std::get<0>(tmp);
            const clade& b = std::get<1>(tmp);

            if (lookup.count(b.data) != 0)
            {
                out.at(i,lookup.at(b.data)) = a + 1;
            }
        }
        cv.clear();
    }

    return out;
}

void mmcutils::median_clade_times(Rcpp::List x, Rcpp::List all_trees, Rcpp::NumericMatrix ts)
{

    Rcpp::List first_tree = x.at(1);
    const int n_tip = rtree(first_tree).n_tip;
    
    const clade_table x_ct(x);
    std::unordered_map<std::vector<bool>, std::vector<double>> time_map;

    for (const clade& cl : x_ct.clades)
    {
        time_map.insert({cl.data, std::vector<double>()});
    }

    const int n_trees = all_trees.length();
    std::vector<std::tuple<int, clade>> cv;    
    cv.reserve(2*n_tip - 1);

    for (int i = 0; i < n_trees; i++)
    {
        Rcpp::List tre_list = all_trees.at(i);
        rtree tr(tre_list);
        find_clades(n_tip, cv, tr);
        for (int j = 0; j < cv.size(); j++)
        {
            const auto& tmp = cv.at(j);
            const int a = std::get<0>(tmp);
            const clade& b = std::get<1>(tmp);

            if (time_map.count(b.data) != 0)
            {
                time_map.at(b.data).push_back(ts.at(i,a));
            }
        }
        cv.clear();
    }

    const auto median = [](auto& v)
    {
        int sz = v.size();
        const int offs = sz/2;

        double out;
        if(sz % 2 == 0)
        {
            std::nth_element(v.begin(), v.begin() + offs, v.end());
            out = v.at(offs) * 0.5;
            std::nth_element(v.begin(), v.begin() + offs - 1, v.end());
            out += v.at(offs-1) * 0.5;
        } else 
        {
            std::nth_element(v.begin(), v.begin() + offs, v.end());
            out = v.at(offs);
        }
        return out;
    };

    std::unordered_map<std::vector<bool>, double> med_map;
    for (auto& it : time_map)
    {
        med_map.insert({it.first, median(it.second)});
    }

    std::vector<double> node_times;

    for (int i = 0; i < x.size(); i++)
    {
        Rcpp::List tre_list = x.at(i);
        rtree tr(tre_list);
        find_clades(n_tip, cv, tr);

        node_times = std::vector<double>(tr.n_node, 0.0);
        for (int j = 0; j < cv.size(); j++)
        {
            const auto& tmp = cv.at(j);
            const int a = std::get<0>(tmp);
            const clade& b = std::get<1>(tmp);

            if (med_map.count(b.data) != 0)
            {
                node_times.at(a) = med_map.at(b.data);
            } else
            {
                throw std::logic_error("Clade not found!");
            }
        }
        cv.clear();
        Rcpp::IntegerMatrix edge = tre_list["edge"]; 
        Rcpp::NumericVector edge_len = tre_list["edge.length"];
        for (int j = 0; j < edge.nrow(); j++)
        {
            double pa_time = node_times.at(edge.at(j,0)-1);
            double ch_time = node_times.at(edge.at(j,1)-1);
            edge_len.at(j) = pa_time - ch_time;
        }            
    }   
}



