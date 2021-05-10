//
// Created by malmiski on 5/2/21.
//

#ifndef BIPARTITE_MATCH_BIPARTITE_MATCH_H
#define BIPARTITE_MATCH_BIPARTITE_MATCH_H
#include <map>
#include <unordered_set>
#include <functional>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Hungarian.h"
using namespace std;
namespace bipartite {
    using point = tuple<double, double>;
    using match = tuple<point, point>;
    using distance_function = function<double(point, point)>;
    using path = vector<point>;
    using edge = tuple<int, int, int>;
    int scaled(double distance, double distance_scale);
    double _log(double a, double base);
    template<typename T>
    unordered_set<T> intersect(unordered_set<T> a, unordered_set<T> b);

    template<typename T>
    unordered_set<T> symmetric_difference(unordered_set<T> &a, unordered_set<T> &b);
    tuple<double, unordered_set<match>> hungarian_algorithm(unordered_set<point> &A, unordered_set<point> &B, const distance_function &dist, double distance_scale);
    unordered_set<match> augment(unordered_set<match> &matches, path &aug_path, unordered_set<point> &A, map<point, int> &point_to_index, vector<int> &weights);
    void PATH_DFS_HELPER(point &start, path &_path, vector<path> &paths, map<int, unordered_set<edge>> &bipartite_graph, vector<point> &A,  unordered_set<point> &free_A, map<point, point> &match_map, map<point, int> &pair_to_index, unordered_set<point> &visited, unordered_set<edge> &eligible_edges);
    void Path_DFS(map<int, unordered_set<edge>> &bipartite_graph, vector<path> &paths, vector<point> &A, vector<point> &B, unordered_set<edge> &eligible_edges, unordered_set<match> &matches, map<point, int> &point_to_index);
    unordered_set<edge> eligible(map<int, unordered_set<edge>> &bipartite_graph, vector<point> &A, vector<point> &B, vector<int> &weights, unordered_set<match> &matches);
    void ADJUST_WEIGHTS_DFS(int _point, unordered_set<int> &visited, vector<int> &weights, unordered_set<int> &S, unordered_set<int> &A,  map<int, unordered_set<edge>> &bipartite_graph, unordered_set<edge> &eligible_edges);
    void adjust_weights(vector<int> &weights, map<int, unordered_set<edge>> &bipartite_graph, unordered_set<point> &free_B, unordered_set<point> &free_A, unordered_set<edge> &eligible_edges, map<point, int> &point_to_index);
    map<int, unordered_set<edge>> create_graph(vector<point> &A, vector<point> &B, distance_function &dist, double distance_scale);
    tuple<unordered_set<point>, unordered_set<point>, unordered_set<match>> Phase_Match(unordered_set<point> &A, unordered_set<point> &B, distance_function &dist, double distance_scale, int stages);
    tuple<double, unordered_set<match>> bipartite_match_edge_cost_scaling(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta, double omega);
    tuple<double, unordered_set<match>> bipartite_match_edge_cost_scaling(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta);
    tuple<double, unordered_set<match>> __bipartite_match_edge_cost_scaling(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta, double omega);
    double matrix_entry(point &a, point &b, vector<vector<double>> matrix);
    tuple<unordered_set<point>, unordered_set<point>, distance_function> convert_matrix_to_points(vector<vector<double>> &cost_matrix);






    }
#include <tuple>
// function has to live in the std namespace
// so that it is picked up by argument-dependent name lookup (ADL).
namespace std{
    namespace
    {

        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        //     and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        //     https://stackoverflow.com/questions/4948780

        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
            static void apply(size_t& seed, Tuple const& tuple)
            {
                HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
            }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
            static void apply(size_t& seed, Tuple const& tuple)
            {
                hash_combine(seed, get<0>(tuple));
            }
        };
    }

    template <typename ... TT>
    struct hash<std::tuple<TT...>>
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
        }

    };
}

#endif //BIPARTITE_MATCH_BIPARTITE_MATCH_H
