#include <map>
#include <set>
#include <functional>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include "Hungarian.h"
#include "bipartite_match.h"
using namespace std;
namespace bipartite{
  int scaled(double distance, double distance_scale){
      return int(ceil(distance*distance_scale));
  }
  double _log(double a, double base){
    return log(a)/log(base);
  }

 template<typename T>
 set<T> intersect(set<T> a, set<T> b){
   set<T> int_set;
   set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(int_set, int_set.begin()));
   return int_set;
 }

    template<typename T>
    set<T> difference(set<T> a, set<T> b) {
        set<T> diff_set;
        set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(diff_set, diff_set.begin()));
        return diff_set;
    }
    template<typename T>
    set<T> symmetric_difference(set<T> a, set<T> b) {
        set<T> diff_set;
        set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(diff_set, diff_set.begin()));
        return diff_set;
    }

    tuple<double, set<match>> hungarian_algorithm(set<point> A, set<point> B, const distance_function &dist, double distance_scale){
    // Create a matrix of each (a,b) with the distance dist(a,b) as the value
    vector<vector<double>> cost_matrix;
    vector<point> vector_A(A.size());
    vector<point> vector_B(B.size());
    copy(A.begin(), A.end(), vector_A.begin());
    copy(B.begin(), B.end(), vector_B.begin());
    for(point _a : A) {
        vector<double> cost(A.size());
        for(int i = 0; i<vector_B.size(); i++){
            cost[i] = scaled(dist(_a, vector_B[i]), distance_scale);
        }
        cost_matrix.push_back(cost);
    }
    vector<int> assignment;
    set<match> matches;
    auto hungarianAlgorithm = HungarianAlgorithm();
    double cost = hungarianAlgorithm.Solve(cost_matrix, assignment);
    for(int i = 0; i<assignment.size(); i++) {
        auto m = make_tuple(vector_A[i], vector_B[assignment[i]]);
        matches.insert(m);
    }
    return make_tuple(cost, matches);
    }



  set<match> augment(set<match> &matches, path &aug_path, set<point> &A, map<point, int> &point_to_index, vector<int> &weights){
    set<match> aug_path_set;
    for(int i = 0; i < aug_path.size() - 1; i++){
      aug_path_set.insert(make_tuple(aug_path[i], aug_path[i+1]));
      int index = point_to_index[aug_path[i]];
      if(index >= A.size()){
        // Decrement the weights of vertices in B that are
        // in the augmenting path to maintain 1-feasible matching
        weights[index] -= 1;
      }
    }
    set<match> new_match;
    new_match = symmetric_difference<match>(matches, aug_path_set);
    return new_match;
  }

  // NOte might need to alter how we are counting vertices as visited, in case we need to determine the
  void PATH_DFS_HELPER(point &start, path &_path, vector<path> &paths, map<int, set<edge>> &bipartite_graph, vector<point> &A,  set<point> &free_A, map<point, point> &match_map, map<point, int> &pair_to_index, set<point> &visited, set<edge> &eligible_edges){
    // Mark this vertex
    visited.insert(start);
    // Next scan all adjacent eligible edges for an unmarked vertex
    int point_index =  pair_to_index[start];
    auto edges = bipartite_graph[point_index];
    set<edge> eligible_current = intersect<edge>(edges, eligible_edges);
    for(auto edge_adjacent : eligible_current){
      auto to_vertex = A[get<1>(edge_adjacent)];
      // If we have visited this vertex already continue
      if(visited.count(to_vertex) == 1){
        continue;
      }
      visited.insert(to_vertex);
      // If this vertex is free, then we add the vertex to the path and finish
      if(free_A.count(to_vertex) == 1){
        _path.push_back(to_vertex);
        auto path_copy = _path;
        _path.clear();
        paths.push_back(path_copy);
        return;
      }

      // Otherwise, this vertex is already matched so let's add the match
      _path.push_back(to_vertex);
      _path.push_back(match_map[to_vertex]);
      PATH_DFS_HELPER(match_map[to_vertex], _path, paths, bipartite_graph, A, free_A, match_map, pair_to_index, visited, eligible_edges);
  }
  // If we didn't get a path from this vertex then remove it and its predecessor
  _path.pop_back();
  if(!_path.empty()) {
      _path.pop_back();
  }
}

  void Path_DFS(map<int, set<edge>> &bipartite_graph, vector<path> &paths, vector<point> &A, vector<point> &B, set<edge> &eligible_edges, set<match> &matches, map<point, int> &point_to_index){
    // Generate set of free vertices of A and B
    set<point> free_A;
    set<point> free_B;
    map<point, point> match_map;
    copy(A.begin(), A.end(), inserter(free_A, free_A.begin()));
    copy(B.begin(), B.end(), inserter(free_B, free_B.begin()));
    for(match m : matches){
      point from = get<0>(m);
      point to = get<1>(m);
      if(point_to_index[from] < A.size()){
        free_A.erase(from);
        free_B.erase(to);
        match_map[from] = to;
      }else{
        free_A.erase(to);
        free_B.erase(from);
        match_map[to] = from;
      }
    }
    set<point> visited;

    while(!free_B.empty()){
      // Start at a free vertex of B, mark this vertex
      point start = *free_B.begin();
      free_B.erase(free_B.begin());
      if(visited.count(start) == 1){
        continue;
      }
      path _path;
      _path.push_back(start);
      PATH_DFS_HELPER(start, _path, paths, bipartite_graph, A, free_A, match_map, point_to_index, visited, eligible_edges);
    }
    // If the unmarked vertex is not free, add the vertex, mark it, and then choose an unmarked vertex
    // from the adjacent eligible edges to find the next starting vertex
    // If there is no such vertex, remove the two recently added vertices and continue
  }
  set<edge> eligible(map<int, set<edge>> &bipartite_graph, vector<point> &A, vector<point> &B, vector<int> &weights, set<match> &matches){
    set<edge> eligible_edges;
    for(int i = A.size(); i<A.size() + B.size(); i++){
      set<edge> _edges = bipartite_graph[i];
      for(edge e : _edges) {
          int from = get<0>(e);
          int to = get<1>(e);
          int dist = get<2>(e);
          bool matched = matches.count(make_tuple(B[from - A.size()], A[to]));
          if (weights[from] + weights[to] == dist && matched) {
              eligible_edges.insert(e);
          } else if (weights[from] + weights[to] == dist + 1 && !matched) {
              eligible_edges.insert(e);
          }
      }
    }
    return eligible_edges;
  }
  void ADJUST_WEIGHTS_DFS(int _point, set<int> &visited, vector<int> &weights, set<int> &S, set<int> &A,  map<int, set<edge>> &bipartite_graph, set<edge> &eligible_edges){
    set<edge> horizon = intersect<edge>(bipartite_graph[_point], eligible_edges);
    for(edge _edge : horizon){
      int to = get<1>(_edge);
      int weight = get<2>(_edge);
      bool matched =  weights[to] + weights[_point] == weight;
      // We will consider the bipartite_graph a directed graph
      // here, where a -> b if (a,b) \in matches, otherwise b -> a
      if((to < A.size() && matched) || (to >= A.size() && !matched)){
        continue;
      }
      if(visited.count(to) == 1){
        continue;
      }
      visited.insert(to);
      if(to < A.size()){
        A.insert(to);
      }else{
        S.insert(to);
      }
      ADJUST_WEIGHTS_DFS(to, visited, weights, S, A, bipartite_graph, eligible_edges);
    }
  }
  void adjust_weights(vector<int> &weights, map<int, set<edge>> &bipartite_graph, set<point> &free_B, set<point> &free_A,
                      set<edge> &eligible_edges, map<point, int> &point_to_index) {
    set<int> B_reachable;
    set<int> A_reachable;
    for(point _point : free_B){
      set<int> B_temp;
      set<int> A_temp;
      set<int> visited;
      // Add the current point to ensure it is counted
      B_reachable.insert(point_to_index[_point]);
      ADJUST_WEIGHTS_DFS(point_to_index[_point], visited, weights, B_temp, A_temp, bipartite_graph, eligible_edges);
      B_reachable.insert(B_temp.begin(), B_temp.end());
      A_reachable.insert(A_temp.begin(), A_temp.end());
    }
    for(int _b : B_reachable){
      weights[_b] += 1;
    }
    for(int _a : A_reachable){
      weights[_a] -= 1;
    }
  }

  map<int, set<edge>> create_graph(vector<point> &A, vector<point> &B, distance_function &dist, double distance_scale){
    map<int, set<edge>> bipartite_graph;
    for(unsigned long i = 0; i<A.size(); i++){
      set<edge> edges;
      for(auto j = 0; j<B.size(); j++){
        // Create an edge with the index of
        edges.insert(make_tuple(i, A.size() + j, scaled(dist(A[i], B[j]), distance_scale)));

        // Add the reverse edge to the other vertex
        if(bipartite_graph.count(int(A.size() + j)) == 0){
          set<edge> new_edges;
          bipartite_graph[A.size() + j] = new_edges;
        }
        bipartite_graph[A.size() + j].insert(make_tuple(A.size() + j, i, scaled(dist(A[i], B[j]), distance_scale)));
      }
      bipartite_graph[int(i)] = edges;
    }
    return bipartite_graph;
  }

  tuple<set<point>, set<point>, set<match>> Phase_Match(set<point> &A, set<point> &B, distance_function &dist, double distance_scale, int stages){

    vector<point> vector_A(A.size());
    vector<point> vector_B(B.size());

    copy(A.begin(), A.end(), vector_A.begin());
    copy(B.begin(), B.end(), vector_B.begin());
    vector<int> weights(A.size() + B.size());

    // Generate dictionary to convert points to their index
    map<point, int> point_to_index;
    for(int i = 0; i < A.size(); i++){
      point_to_index[vector_A[i]] = i;
    }

    for(int i = 0; i < B.size(); i++){
      point_to_index[vector_B[i]] = int(A.size() + i);
    }
    // At first all vertices are free and matches are empty, as we have not begun matching yet
    set<match> matches;
    set<point> free_A;
    set<point> free_B;
    // Create the bipartite graph based on the point sets, using array indices
    map<int, set<edge>> bipartite_graph = create_graph(vector_A, vector_B, dist, distance_scale);

    for(int k = 0; k<stages; k++){
      // Generate set of eligible edges
      set<edge> eligible_edges = eligible(bipartite_graph, vector_A, vector_B, weights, matches);
      vector<path> paths;
      Path_DFS(bipartite_graph, paths, vector_A, vector_B, eligible_edges, matches, point_to_index);
      for(path aug_path : paths){
        matches = augment(matches, aug_path, A, point_to_index, weights);
      }
      copy(A.begin(), A.end(), inserter(free_A, free_A.begin()));
      copy(B.begin(), B.end(), inserter(free_B, free_B.begin()));
      for(auto [from, to] : matches){
          free_B.erase(from);
          free_A.erase(to);
      }
        adjust_weights(weights, bipartite_graph, free_B, free_A, eligible_edges, point_to_index);
      // If the matches are a complete matching exit
      if(matches.size() == A.size()){
        break;
      }
    }
    return make_tuple(free_A, free_B, matches);
  }

  tuple<double, set<match>> bipartite_match_edge_cost_scaling(const set<point> &A, const set<point> &B,  distance_function &dist, double delta, double omega){
    int match_size = 0;
    set<match> matches;
    double cost = 0;
    int n = int(A.size());
    int n_two_thirds = int(pow(n, 2.0/3.0));
    set<point> A_i = A;
    set<point> B_i = B;
    double epsilon = 1/(2*_log(1/delta, 3));
    double distance_scale = 2*n/(epsilon * omega);
    int i = 0;
    while(match_size < n){
        int n_i = int(A_i.size());
        int stages = int(ceil((30/epsilon)*pow(n, pow(3, i)*delta)));
        if(n_i < n_two_thirds){
            auto [hungarian_cost, hungarian_matches] = hungarian_algorithm(A_i, B_i, dist, distance_scale);
            matches.insert(hungarian_matches.begin(), hungarian_matches.end());
            cost += hungarian_cost;
            return make_tuple(cost, matches);
        }else{
            auto [A_i_new, B_i_new, new_matches] = Phase_Match(A_i, B_i, dist, distance_scale, stages);
            matches.insert(new_matches.begin(), new_matches.end());
            distance_scale *= 1/(2*(1 + epsilon)*(1+epsilon)*pow(n, (pow(3,i)*delta) - 1));
            match_size += int(new_matches.size());
            for(auto new_match : new_matches){
                cost += dist(get<1>(new_match), get<0>(new_match));
            }
            A_i = A_i_new;
            B_i = B_i_new;
        }
        i += 1;
    }
    return make_tuple(cost, matches);
}
    tuple<double, set<match>> bipartite_match_edge_cost_scaling(const set<point> &A, const set<point> &B,  distance_function &dist, double delta){
        const double n_approximate_match = 35484;// Add a call to an n-approximate matching for A and B that runs in O(n^2) time
        double min_cost = INFINITY;
        set<match> min_match;
        for(int i = 0; i<ceil(_log(A.size(), 2)); i++){
            auto [new_cost, matches] = bipartite_match_edge_cost_scaling(A, B, dist, delta, n_approximate_match/pow(2, i));
            if(min_cost > new_cost){
                min_cost = new_cost;
                min_match = matches;
            }
        }
        return make_tuple(min_cost, min_match);
    }

    double matrix_entry(point &a, point &b, vector<vector<double>> matrix){
      return matrix[get<0>(a)-1][get<1>(b)-1];
    }
    tuple<set<point>, set<point>, distance_function> convert_matrix_to_points(vector<vector<double>> &cost_matrix){
      // Only for matrices of size n x n
      set<point> A;
      set<point> B;

      for(int i = 0; i<cost_matrix.size(); i++){
          A.insert(make_tuple(i+1, 0));
          B.insert(make_tuple(0, i+1));
      }
      auto dist = [cost_matrix](point a, point b) -> double{return matrix_entry(a, b, cost_matrix);};
      return make_tuple(A, B, dist);
  }
}
