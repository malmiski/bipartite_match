#include "bipartite_match.cpp"
int main() {
    std::cout << "Hello, World!" << std::endl;
    vector<vector<double> > costMatrix = {{10, 19, 8, 15, 0},
                                          {10, 18, 7, 17, 0},
                                          {13, 16, 9, 14, 0},
                                          {12, 19, 8, 18, 0}};

    HungarianAlgorithm HungAlgo;
    vector<int> assignment;

    double cost = HungAlgo.Solve(costMatrix, assignment);
    auto [A, B, dist] = bipartite::convert_matrix_to_points(costMatrix);
    bipartite::bipartite_match_edge_cost_scaling(A, B, dist, 0.01);
    for (unsigned int x = 0; x < costMatrix.size(); x++)
        std::cout << x << "," << assignment[x] << "\t";

    std::cout << "\ncost: " << cost << std::endl;

    return 0;
}