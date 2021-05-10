#include <string>
#include <unordered_set>
#include <cmath>
#include <tuple>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "bipartite_match.h"
#include <chrono>
#include "bipartite_match_ANN.h"
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
    srand(100); // Seed random to known constant
    vector<double> costs(10);
    for(int i = 1; i<10; i++){
        char buff[100];
        snprintf(buff, sizeof(buff), "../tests/bipartite_set_%d.txt", i);
        auto [A_, B_, dist_] = load(buff);
        auto current = chrono::system_clock::now();
        auto [cost, hungarianMatches] = bipartite::hungarian_algorithm(A_, B_, dist_, 1);//HungAlgo.Solve(costMatrix, assignment);
        auto end = chrono::system_clock::now();
        chrono::duration<double> duration = end - current;
        cout << "Bipartite set " << i << ":" << endl;
        cout << "Time: " << duration.count() << " seconds\n";
        cout << "cost: " << cost << endl << endl;
        costs.push_back(cost);
    }
    auto deltas = {.08, .15, .3};
    for(auto delta: deltas)
    for(int i = 1; i<=10; i++){
        char buff[100];
        snprintf(buff, sizeof(buff), "../tests/bipartite_set_%d.txt", i);
        auto [A_, B_, dist_] = load(buff);
        auto current = chrono::system_clock::now();
        auto [cost, hungarianMatches] = bipartite::__bipartite_match_edge_cost_scaling(A_, B_, dist_, delta, costs[i]*3/4);
        auto end = chrono::system_clock::now();
        chrono::duration<double> duration = end - current;
        cout << "(\\delta = " << delta << ") Bipartite set " << i << ":" << endl;
        cout << "Time: " << duration.count() << " seconds\n";
        cout << "cost: " << cost << endl << endl;
    }
    deltas = {.15, .3};
    auto varepsilons = {.75, .25};
    for(auto delta: deltas)
    for(auto varepsilon : varepsilons)
    for(int i = 1; i<=10; i++){
        char buff[100];
        snprintf(buff, sizeof(buff), "../tests/bipartite_set_%d.txt", i);
        auto [A_, B_, dist_] = load(buff);
        auto current = chrono::system_clock::now();
        auto [cost, hungarianMatches] = __bipartite_match_edge_cost_scaling_ANN(A_, B_, dist_, delta, costs[i]*3/4, varepsilon);//HungAlgo.Solve(costMatrix, assignment);
        auto end = chrono::system_clock::now();
        chrono::duration<double> duration = end - current;
        cout << "(\\delta = " << delta << ", \\varepsilon = " << varepsilon << ")Bipartite set " << i << ":" << endl;
        cout << "Time: " << duration.count() << " seconds\n";
        cout << "cost: " << cost << endl << endl;
        costs.push_back(cost);
    }

    // In the paper they use an n^2 approximation algorithm to get an n-approximate matching to use in this
    // algorithm, however, I couldn't find what algorithm they used so instead I will be using the hungarian
    // algorithms' cost multiplied by a random number between 1 and n to get the same effect
//    auto random_cost = (1 + ((double(rand()) / RAND_MAX) * double(A_.size()-1))) * cost;
//    auto [approximate_cost, approximate_matches] = __bipartite_match_edge_cost_scaling_ANN(A_, B_, dist_, 0.14, cost*2/3,.6);
//    std::cout << "\nApproximation cost: " << approximate_cost << std::endl;
//    end = std::chrono::system_clock::now();
//    duration = end - current;
//    std::cout << "Took: " << duration.count() << " seconds\n";
//    cout << endl;
    return 0;
}