//
// Created by malmiski on 5/2/21.
//

#ifndef BIPARTITE_MATCH_BIPARTITE_MATCH_H
#define BIPARTITE_MATCH_BIPARTITE_MATCH_H
namespace bipartite {
    using point = tuple<double, double>;
    using match = tuple<point, point>;
    using distance_function = function<double(point, point)>;
    using path = vector<point>;
    using edge = tuple<int, int, int>;

    int scaled(double distance, double distance_scale);
    double _log(double a, double base);
    template<typename T>
    set<T> intersect(set<T> a, set<T> b);

    template<typename T>
    set<T> difference(set<T> a, set<T> b);
    template<typename T>
    set<T> symmetric_difference(set<T> a, set<T> b);
    tuple<double, set<match>> hungarian_algorithm(set<point> A, set<point> B, const distance_function &dist, double distance_scale);
    set<match> augment(set<match> &matches, path &aug_path, set<point> &A, map<point, int> point_to_index, vector<int> weights);
    void PATH_DFS_HELPER(point &start, path &_path, vector<path> &paths, map<int, set<edge>> &bipartite_graph, vector<point> &A,  set<point> &free_A, map<point, point> &match_map, map<point, int> &pair_to_index, set<point> &visited, set<edge> &eligible_edges);
    void Path_DFS(map<int, set<edge>> &bipartite_graph, vector<path> &paths, vector<point> &A, vector<point> &B, set<edge> &eligible_edges, set<match> &matches, map<point, int> &point_to_index);
    b






    }
#endif //BIPARTITE_MATCH_BIPARTITE_MATCH_H
