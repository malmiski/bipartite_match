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
    set<match> augment(set<match> &matches, path &aug_path, set<point> &A, map<point, int> &point_to_index, vector<int> &weights);
    void PATH_DFS_HELPER(point &start, path &_path, vector<path> &paths, map<int, set<edge>> &bipartite_graph, vector<point> &A,  set<point> &free_A, map<point, point> &match_map, map<point, int> &pair_to_index, set<point> &visited, set<edge> &eligible_edges);
    void Path_DFS(map<int, set<edge>> &bipartite_graph, vector<path> &paths, vector<point> &A, vector<point> &B, set<edge> &eligible_edges, set<match> &matches, map<point, int> &point_to_index);
    set<edge> eligible(map<int, set<edge>> &bipartite_graph, vector<point> &A, vector<point> &B, vector<int> &weights, set<match> &matches);
    void ADJUST_WEIGHTS_DFS(int _point, set<int> &visited, vector<int> &weights, set<int> &S, set<int> &A,  map<int, set<edge>> &bipartite_graph, set<edge> &eligible_edges);
    void adjust_weights(vector<int> &weights, map<int, set<edge>> &bipartite_graph, set<point> &free_B, set<point> &free_A, set<edge> &eligible_edges, map<point, int> &point_to_index);
    map<int, set<edge>> create_graph(vector<point> &A, vector<point> &B, distance_function &dist, double distance_scale);
    tuple<set<point>, set<point>, set<match>> Phase_Match(set<point> &A, set<point> &B, distance_function &dist, double distance_scale, int stages);
    tuple<double, set<match>> bipartite_match_edge_cost_scaling(const set<point> &A, const set<point> &B,  distance_function &dist, double delta, double omega);
    tuple<double, set<match>> bipartite_match_edge_cost_scaling(const set<point> &A, const set<point> &B,  distance_function &dist, double delta);
    double matrix_entry(point &a, point &b, vector<vector<double>> matrix);
    tuple<set<point>, set<point>, distance_function> convert_matrix_to_points(vector<vector<double>> &cost_matrix);






    }
#endif //BIPARTITE_MATCH_BIPARTITE_MATCH_H
