//
// Created by malmiski on 5/8/21.
//

#include "bipartite_match_ANN.h"

namespace bipartite{
    template<typename T>
    unordered_set<T> symmetric_difference_ANN(unordered_set<T> &a, unordered_set<T> &b) {
        unordered_set<T> diff_set(a.size() + b.size());
        while(!a.empty()){
            auto a_s = *a.begin();
            a.erase(a_s);
            if(b.count(a_s) == 1){
                a.erase(a_s);
                b.erase(a_s);
            }else{
                diff_set.insert(a_s);
            }
        }
        diff_set.insert(b.begin(), b.end());
        return diff_set;
    }

    unordered_set<match> augment_ANN(unordered_set<match> &matches, path &aug_path, unordered_set<point> &A, map<point, int> &point_to_index, vector<int> &weights, distance_function &dist, double distance_scale){
        unordered_set<match> aug_path_set;
        for(int i = 0; i < aug_path.size() - 1; i++){
            int index = point_to_index[aug_path[i]];
            if(index >= A.size()) aug_path_set.insert(make_tuple(aug_path[i], aug_path[i+1]));
            else aug_path_set.insert(make_tuple(aug_path[i+1], aug_path[i]));
            if(index >= A.size()){
                // Decrement the weights of vertices in B that are
                // in the augmenting path to maintain epsilon-feasible matching
                weights[index] = scaled(dist(aug_path[i], aug_path[i+1]), distance_scale) - weights[point_to_index[aug_path[i+1]]];
            }
        }
        unordered_set<match> new_match;
        new_match = symmetric_difference_ANN<match>(matches, aug_path_set);
        return new_match;
    }

    // NOte might need to alter how we are counting vertices as visited, in case we need to determine the
    bool PATH_DFS_HELPER_ANN(point &start, path &_path, vector<path> &paths, map<int, unordered_set<edge>> &bipartite_graph, vector<point> &A,  unordered_set<point> &free_A, map<point, point> &match_map, map<point, int> &pair_to_index, ANN_DS &ann_ds, unordered_set<edge> &eligible_edges){
        // Next scan all adjacent eligible edges for an unmarked vertex
        int point_index =  pair_to_index[start];
        auto edges = bipartite_graph[point_index];
        edge edge_adjacent;

        while(ann_ds.eligible(start, edge_adjacent)){
            eligible_edges.insert(edge_adjacent);
            auto to_vertex = A[get<1>(edge_adjacent)];
            ann_ds.delete_point(get<1>(edge_adjacent));
            // If this vertex is free, then we add the vertex to the path and finish
            if(free_A.count(to_vertex) == 1){
                _path.push_back(to_vertex);
                auto path_copy = _path;
                _path.clear();
                paths.push_back(path_copy);
                return true;
            }

            // Otherwise, this vertex is already matched so let's add the match
            _path.push_back(to_vertex);
            _path.push_back(match_map[to_vertex]);
            if(PATH_DFS_HELPER_ANN(match_map[to_vertex], _path, paths, bipartite_graph, A, free_A, match_map, pair_to_index, ann_ds, eligible_edges)){
                return true;
            }
        }
        // If we didn't get a path from this vertex then remove it and its predecessor
        _path.pop_back();
        if(!_path.empty()) {
            _path.pop_back();
        }
        return false;
    }

    void Path_DFS_ANN(map<int, unordered_set<edge>> &bipartite_graph, vector<path> &paths, vector<point> &A, vector<point> &B, ANN_DS &ann_ds, unordered_set<match> &matches, map<point, int> &point_to_index, unordered_set<edge> &eligible_edges){
        // Generate set of free vertices of A and B
        unordered_set<point> free_A(A.begin(), A.end());
        unordered_set<point> free_B(B.begin(), B.end());
        map<point, point> match_map;
//    copy(A.begin(), A.end(), inserter(free_A, free_A.begin()));
//    copy(B.begin(), B.end(), inserter(free_B, free_B.begin()));
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

        while(!free_B.empty()){
            // Start at a free vertex of B, mark this vertex
            point start = *free_B.begin();
            free_B.erase(free_B.begin());
            path _path;
            _path.push_back(start);
            PATH_DFS_HELPER_ANN(start, _path, paths, bipartite_graph, A, free_A, match_map, point_to_index, ann_ds, eligible_edges);
        }
        // If the unmarked vertex is not free, add the vertex, mark it, and then choose an unmarked vertex
        // from the adjacent eligible edges to find the next starting vertex
        // If there is no such vertex, remove the two recently added vertices and continue
    }

    tuple<unordered_set<point>, unordered_set<point>, unordered_set<match>> Phase_Match_ANN(unordered_set<point> &A, unordered_set<point> &B, distance_function &dist, double distance_scale, int stages, double varepsilon, ANN_DS &ann_ds){

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
        unordered_set<match> matches;
        unordered_set<point> free_A;
        unordered_set<point> free_B;
        // Create the bipartite graph based on the point sets, using array indices
        map<int, unordered_set<edge>> bipartite_graph = create_graph(vector_A, vector_B, dist, distance_scale);
        ann_ds.A = vector_A;
        ann_ds.weights = weights;
        ann_ds.point_to_index = point_to_index;
        ann_ds.d = [distance_scale, dist](auto a, auto b) -> double{return int(ceil(distance_scale*dist(a,b)));};
        free_A = unordered_set<point>(A.begin(), A.end());
        free_B = unordered_set<point>(B.begin(), B.end());
        for(int k = 0; k<stages; k++){
            ann_ds.set_current_weight(k);
            // Generate set of eligible edges
            unordered_set<edge> eligible_edges;
            vector<path> paths;
            Path_DFS_ANN(bipartite_graph, paths, vector_A, vector_B, ann_ds, matches, point_to_index, eligible_edges);
            for(auto &aug_path : paths){
                matches = augment_ANN(matches, aug_path, A, point_to_index, weights, dist, distance_scale);
                for(auto [from, to] : matches){
                    free_B.erase(from);
                    free_A.erase(to);
                }
                auto a = 0;
                auto b = a;
            }
            adjust_weights(weights, bipartite_graph, free_B, free_A, eligible_edges, point_to_index);
            // If the matches are a complete matching exit
            if(matches.size() == A.size()){
                break;
            }
        }
        return make_tuple(free_A, free_B, matches);
    }

    tuple<double, unordered_set<match>> __bipartite_match_edge_cost_scaling_ANN(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta, double omega, double varepsilon){
        int match_size = 0;
        unordered_set<match> matches;
        double cost = 0;
        int n = int(A.size());
        int n_one_thirds = int(pow(n, 1.0 / 3.0));
        unordered_set<point> A_i = A;
        unordered_set<point> B_i = B;
        double epsilon = 1/(2*_log(1/delta, 3));
        double distance_scale = 2*n/(epsilon * omega);
        int i = 0;
        auto ann_ds = ANN_DS(varepsilon, n);
        while(match_size < n){
            int n_i = int(A_i.size());
            int stages = int(ceil((30/epsilon)*pow(n, pow(2, i)*delta)));
            if(n_i < n_one_thirds){
                auto [hungarian_cost, hungarian_matches] = hungarian_algorithm(A_i, B_i, dist, distance_scale);
                matches.insert(hungarian_matches.begin(), hungarian_matches.end());
                cost += hungarian_cost;
                return make_tuple(cost, matches);
            }else{
                auto [A_i_new, B_i_new, new_matches] = Phase_Match_ANN(A_i, B_i, dist, distance_scale, stages, varepsilon, ann_ds);
                matches.insert(new_matches.begin(), new_matches.end());
                distance_scale *= 1/(2*(1 + epsilon)*(1+epsilon)*pow(n, (pow(2,i-1)*delta)));
                match_size += int(new_matches.size());
                for(auto new_match : new_matches){
                    cost += ceil(dist(get<1>(new_match), get<0>(new_match)));
                }
                A_i = A_i_new;
                B_i = B_i_new;
            }
            i += 1;
        }
        return make_tuple(cost, matches);
    }
    tuple<double, unordered_set<match>> bipartite_match_edge_cost_scaling_ANN(const unordered_set<point> &A, const unordered_set<point> &B,  distance_function &dist, double delta, double approximate_cost, double varepsilon){
        const double n_approximate_match = approximate_cost;// Add a call to an n-approximate matching for A and B that runs in O(n^2) time
        double min_cost = INFINITY;
        unordered_set<match> min_match;
        for(int i = 1; i<=int(ceil(_log(double(A.size()), 2)))-1; i++){
            auto [new_cost, matches] = __bipartite_match_edge_cost_scaling_ANN(A, B, dist, delta, n_approximate_match/pow(2, i), varepsilon);
            if(min_cost > new_cost){
                min_cost = new_cost;
                min_match = matches;
            }
        }
        return make_tuple(min_cost, min_match);
    }


    bipartite::ANN_DS::ANN_DS(double varepsilon, int n):varepsilon(varepsilon), n(n){
        build_wspd(int(ceil(5 / varepsilon * n)), varepsilon / 10);
    }

    bool bipartite::ANN_DS::eligible(point &b, edge &_edge) {
        auto weight = this->current_weight;
        auto nn = new ANNidx [1];
        auto dd = new ANNdist[1];
        double p[] = {get<0>(b), get<1>(b)};

        for(auto tree : current_trees){
            if(tree == nullptr)
                continue;
            tree->annkSearch(p, 1, nn, dd, varepsilon/2);
            auto a = tree->thePoints()[nn[0]];
//            delete[] nn;
//            delete[] dd;
            auto a_point = make_tuple(a[0], a[1]);
            int index = point_to_index[a_point];
            auto sum = weights[index] + this->current_weight;
            auto distance = d(a_point, b);
            auto epsilon_distance = ceil(varepsilon*distance);
            if(distance + 1 <= sum && sum <= distance + epsilon_distance){
                _edge = make_tuple(point_to_index[b], index, d(a_point,b));
                return true;
            }
        }
        return false;
    }
    void bipartite::ANN_DS::delete_point(int index) {
        point a = A[index];
        for(int i = 0; i<K_s.size(); i++){
            if(K_s[i].count(a) == 1){
                K_s[i].erase(a);
                delete current_trees[i];
                current_trees[i] = build_ANN(K_s[i]);
            }
        }
    }
    unordered_set<point> bipartite::ANN_DS::get_points(int left, int right) {
        unordered_set<point> K;
        for(int i =0; i<A.size(); i++){
            if(left <= -weights[i] && -weights[i] <= right){
                K.insert(A[i]);
            }
        }
        return K;
    }

    void bipartite::ANN_DS::set_current_weight(int weight){
        this->current_weight = weight;
        current_partitions.clear();
        current_trees.clear();
        K_s.clear();
        for(auto interval : partitions){
            auto y = get<1>(interval);
            auto l = get<0>(y);
            auto r = get<1>(y);
            if(weight <= r && weight >= l){
                auto x = get<0>(interval);
                current_partitions.emplace_back(interval);
                auto K = get_points(get<0>(x), get<1>(x));
                if(K.empty()) {
                    current_trees.emplace_back(nullptr);
                    K_s.emplace_back(unordered_set<point>());
                    continue;
                }
                current_trees.emplace_back(build_ANN(K));
                K_s.emplace_back(K);
            }
        }
    }
    ANNkd_tree* bipartite::ANN_DS::build_ANN(unordered_set<point> &K) {
        auto pts = annAllocPts(int(K.size()), 2);
        int i = 0;
        for(auto k : K){
            pts[i][0] = get<0>(k);
            pts[i][1] = get<1>(k);
            i += 1;
        }
        return new ANNkd_tree(pts, int(K.size()), 2);
    }

    void bipartite::ANN_DS::build_wspd(int N, double rho){
        auto w = bipartite_match::wspd(N, rho);
        partitions = vector<tuple<interval, interval>>(w.intervals.size());
        copy(w.intervals.begin(), w.intervals.end(), partitions.begin());
    }

}