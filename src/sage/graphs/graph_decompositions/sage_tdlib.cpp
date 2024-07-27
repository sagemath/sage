// sage_setup: distribution = sagemath-tdlib

#include <boost/tuple/tuple.hpp>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <treedec/combinations.hpp>
#include <treedec/misc.hpp>

#ifndef TD_STRUCT_VERTEX
#define TD_STRUCT_VERTEX

struct Vertex{
    unsigned int id;
};

#endif

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex> TD_graph_t;

void make_tdlib_graph(TD_graph_t &G, std::vector<unsigned int> &V, std::vector<unsigned int> &E){
    unsigned int max = 0;
    for(unsigned int i = 0; i < V.size(); i++)
        max = (V[i]>max)? V[i] : max;

    std::vector<TD_graph_t::vertex_descriptor> idxMap(max+1);

    for(unsigned int i = 0; i < V.size(); i++){
        idxMap[V[i]] = boost::add_vertex(G);
        G[idxMap[V[i]]].id = V[i];
    }

    if(E.size() != 0){
        for(unsigned int j = 0; j < E.size()-1; j++){
            boost::add_edge(idxMap[E[j]], idxMap[E[j+1]], G);
            j++;
        }
    }
}

void make_sage_decomp(TD_tree_dec_t &T, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    std::map<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int> vertex_map;
    boost::graph_traits<TD_tree_dec_t>::vertex_iterator tIt, tEnd;
    unsigned int id = 0;
    
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        vertex_map.insert(std::pair<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int>(*tIt, id++));
        std::vector<int> bag;
        for(std::set<unsigned int>::iterator sIt = T[*tIt].bag.begin(); sIt != T[*tIt].bag.end(); sIt++)
            bag.push_back((int)*sIt);
        V_T.push_back(bag);
    }
    
    boost::graph_traits<TD_tree_dec_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(T); eIt != eEnd; eIt++){
        std::map<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int>::iterator v, w;
        v = vertex_map.find(boost::source(*eIt, T));
        w = vertex_map.find(boost::target(*eIt, T));
        E_T.push_back(v->second);
        E_T.push_back(w->second);
    }
}


/* EXACT TREE DECOMPOSITIONS */

int sage_exact_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::exact_decomposition_cutset(G, T, lb);

    treedec::make_small(T);

    make_sage_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

