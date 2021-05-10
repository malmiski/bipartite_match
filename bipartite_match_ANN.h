//
// Created by malmiski on 5/8/21.
//

#ifndef BIPARTITE_MATCH_BIPARTITE_MATCH_ANN_H
#define BIPARTITE_MATCH_BIPARTITE_MATCH_ANN_H
#include <map>
#include <unordered_set>
#include <functional>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Hungarian.h"
#include "include/ANN/ANN.h"
#include "wspd.h"
#include "bipartite_match.h"

using namespace std;
using namespace bipartite_match;
namespace bipartite {
    class ANN_DS{
        int n;
        double varepsilon;
        double current_weight;
        vector<tuple<interval, interval>> partitions;
        vector<tuple<interval, interval>> current_partitions;
        vector<unordered_set<point>> K_s;
        vector<ANNkd_tree*> current_trees;
        static ANNkd_tree* build_ANN(unordered_set<point> &K);
        void build_wspd(int N, double rho);
        unordered_set<point> get_points(int left, int right);
    public:
        distance_function d;
        vector<int> weights;
        vector<point> A;
        map<point, int> point_to_index;
        ANN_DS(double varepsilon, int n);
        bool eligible(point &b, edge &_return_edge);
        void set_current_weight(int weight);
        void delete_point(int a);
    };
    unordered_set<match> augment_ANN(unordered_set<match> &matches, path &aug_path, unordered_set<point> &A, map<point, int> &point_to_index, vector<int> &weights, distance_function &dist, double distance_scale);
    bool PATH_DFS_HELPER_ANN(point &start, path &_path, vector<path> &paths, map<int, unordered_set<edge>> &bipartite_graph, vector<point> &A,  unordered_set<point> &free_A, map<point, point> &match_map, map<point, int> &pair_to_index, ANN_DS &ann_ds, unordered_set<edge> &eligible_edges);
    void Path_DFS_ANN(map<int, unordered_set<edge>> &bipartite_graph, vector<path> &paths, vector<point> &A, vector<point> &B, ANN_DS &ann_ds, unordered_set<match> &matches, map<point, int> &point_to_index, unordered_set<edge> &eligible_edges);
    tuple<unordered_set<point>, unordered_set<point>, unordered_set<match>> Phase_Match_ANN(unordered_set<point> &A, unordered_set<point> &B, distance_function &dist, double distance_scale, int stages, double varepsilon);
    tuple<double, unordered_set<match>> __bipartite_match_edge_cost_scaling_ANN(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta, double omega, double varepsilon);
    tuple<double, unordered_set<match>> bipartite_match_edge_cost_scaling_ANN(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta, double approximatecost, double varepsilon);
}

#endif //BIPARTITE_MATCH_BIPARTITE_MATCH_ANN_H
