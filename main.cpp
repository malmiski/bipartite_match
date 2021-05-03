#include <string>
#include <cmath>
#include <tuple>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "bipartite_match.cpp"
using namespace std;
using namespace bipartite;
set<point> process(string line){
    auto index = 0;
    set<point> points;
    int n = line.length();
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
            strings.push_back("");
            while(j < line.length() && k != ')'){
                k = line.at(j);
                if(k == ','){
                    strings.push_back("");
                    j += 1;
                    continue;
                }
                strings.back() += k;
                j += 1;
            }
            index = j;
            double second = stod(*strings.begin(), NULL);
            double first = stod(strings.back(), NULL);
            points.insert(make_tuple(first, second));
        }
        index++;
    }
    return points;
}
tuple<set<point>, set<point>, distance_function> load(string filename){
    set<point> A;
    set<point> B;
    ifstream afile(filename);
    string point_string;
    getline (afile, point_string);
    A = process(point_string);
    getline (afile, point_string);
    B = process(point_string);
    afile.close();
    auto dist = [](auto a, auto b) -> double {return sqrt(pow(get<0>(a) - get<0>(b), 2) + pow(get<1>(a) - get<1>(b), 2));};
    return make_tuple(A, B, dist);
}
int main() {
    std::cout << "Hello, World!" << std::endl;
    srand(100); // Seed random to known constant
    auto [A_, B_, dist_] = load("../tests/bipartite_set_1.txt");
    vector<vector<double> > costMatrix = {{10, 19, 8, 15, 0},
                                          {10, 18, 7, 17, 0},
                                          {13, 16, 9, 14, 0},
                                          {12, 19, 8, 18, 0}};

    HungarianAlgorithm HungAlgo;
    vector<int> assignment;

    auto [cost, hungoMatches] = bipartite::hungarian_algorithm(A_, B_, dist_, 1);//HungAlgo.Solve(costMatrix, assignment);
//    auto [A, B, dist] = bipartite::convert_matrix_to_points(costMatrix);
    std::cout << "Hello, World!" << std::endl;
    std::cout << "\ncost: " << cost << std::endl;
    // In the paper they use an n^2 approximation algorithm to get an n-approximate matching to use in this
    // algorithm, however, I couldn't find what algorithm they used so instead I will be using the hungarian
    // algorithms' cost multiplied by a random number between 1 and n to get the same effect
    auto random_cost = (1 + ((double(rand()) / RAND_MAX) * (A_.size()-1))) * cost;
    auto [acost, aMatches] = bipartite::bipartite_match_edge_cost_scaling(A_, B_, dist_, 0.5, random_cost);
//    for (unsigned int x = 0; x < costMatrix.size(); x++)
//        std::cout << x << "," << assignment[x] << "\t";
//
    std::cout << "\nApproximation cost: " << acost << std::endl;

    return 0;
}