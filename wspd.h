//
// Created by malmiski on 5/8/21.
//

#ifndef BIPARTITE_MATCH_WSPD_H
#define BIPARTITE_MATCH_WSPD_H
#include <unordered_set>
#include <tuple>
#include "bipartite_match.h"
using namespace std;
using interval = tuple<int, int>;
namespace bipartite_match {
    class wspd {
    private:
        class node {
        public:
            int l;
            int r;
            node *left, *right;
            node(int left, int right) : l(left), r(right) {}
//            ~node(){if(left) delete left; if(right) delete right;}
            inline int length() {return this->r - this->l;}
        };
        node* construct(int l, int r);
        void __wspd(node*, node*);
        node* split_tree;
        double rho;
        double N;
    public:
        unordered_set<tuple<interval, interval>> intervals;
        wspd(int N, double rho);
        ~wspd(){delete split_tree;}
    };

}
#endif //BIPARTITE_MATCH_WSPD_H
