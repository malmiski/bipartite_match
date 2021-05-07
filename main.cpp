#include <string>
#include <unordered_set>
#include <cmath>
#include <tuple>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "bipartite_match.cpp"
#include <chrono>
using namespace std;
using namespace bipartite;
unordered_set<point> process(string line){
    auto index = 0;
    unordered_set<point> points;
    int n = int(line.length());
    while(index < n){
        auto c = line.at(index);
        if(c == '[' || c == ']'){
            index++;
            continue;
        }
        if(c == '('){
            auto j = index + 1;
            auto k = line.at(j);
            vector<string> strings;
            strings.emplace_back("");
            while(j < line.length() && k != ')'){
                k = line.at(j);
                if(k == ','){
                    strings.emplace_back("");
                    j += 1;
                    continue;
                }
                strings.back() += k;
                j += 1;
            }
            index = j;
            double second = stod(*strings.begin(), nullptr);
            double first = stod(strings.back(), nullptr);
            points.insert(make_tuple(first, second));
        }
        index++;
    }
    return points;
}
tuple<unordered_set<point>, unordered_set<point>, distance_function> load(const string &filename){
    unordered_set<point> A;
    unordered_set<point> B;
    ifstream input_file(filename);
    string point_string;
    getline (input_file, point_string);
    A = process(point_string);
    getline (input_file, point_string);
    B = process(point_string);
    input_file.close();
    auto dist = [](auto a, auto b) -> double {return sqrt(pow(get<0>(a) - get<0>(b), 2) + pow(get<1>(a) - get<1>(b), 2));};
    return make_tuple(A, B, dist);
}
int main() {
    std::cout << "Hello, World!" << std::endl;
    auto current = std::chrono::system_clock::now();
    srand(100); // Seed random to known constant
    auto [A_, B_, dist_] = load("../tests/bipartite_set_2.txt");
//    vector<vector<double> > costMatrix = {{10, 19, 8, 15, 0},
//                                          {10, 18, 7, 17, 0},
//                                          {13, 16, 9, 14, 0},
//                                          {12, 19, 8, 18, 0}};

//    HungarianAlgorithm HungAlgo;
//    vector<int> assignment;

    auto [cost, hungarianMatches] = bipartite::hungarian_algorithm(A_, B_, dist_, 1);//HungAlgo.Solve(costMatrix, assignment);
//    auto [A, B, dist] = bipartite::convert_matrix_to_points(costMatrix);
    std::cout << "Hello, World!" << std::endl;
    std::cout << "\ncost: " << cost << std::endl;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end - current;
    std::cout << "Took: " << duration.count() << " seconds\n";
    current = end;
    // In the paper they use an n^2 approximation algorithm to get an n-approximate matching to use in this
    // algorithm, however, I couldn't find what algorithm they used so instead I will be using the hungarian
    // algorithms' cost multiplied by a random number between 1 and n to get the same effect
    auto random_cost = (1 + ((double(rand()) / RAND_MAX) * double(A_.size()-1))) * cost;
    auto [approximate_cost, approximate_matches] = bipartite::__bipartite_match_edge_cost_scaling(A_, B_, dist_, 0.14, cost*2/3.0);
    std::cout << "\nApproximation cost: " << approximate_cost << std::endl;
    end = std::chrono::system_clock::now();
    duration = end - current;
    std::cout << "Took: " << duration.count() << " seconds\n";
    return 0;
}