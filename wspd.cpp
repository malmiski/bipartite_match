//
// Created by malmiski on 5/8/21.
//

#include "wspd.h"
using namespace bipartite_match;
wspd::node* wspd::construct(int l, int r){
    if(l == r) {
        auto n = new node(l, r);
        return n;
    }
    auto root = new node(l, r);
    auto middle = l + (r - l)/2;
    root->left = construct(l, middle);
    root->right = construct(middle + 1, r);
    return root;
}

wspd::wspd(int N, double rho) {
    this->split_tree = construct(0, N);
    this->N = N;
    this->rho = rho;
    __wspd(this->split_tree, this->split_tree);
}

void wspd::__wspd(wspd::node *l, wspd::node *r) {
    if(l->length() < r->length()){
        auto tmp = l;
        l = r;
        r = tmp;
    }
    auto dist = 0;
    if(l->r < r->l){
        dist = r->l - l->r;
    }else{
        dist = l->l - r->r;
    }
    if(l->length() <= this->rho*dist){
        if(l->r > r->l){
            auto tmp = l;
            l = r;
            r = tmp;
        }
        if(l->r == r->r and l->l == r->l){
            return;
        }
        this->intervals.insert(make_tuple(make_tuple(l->l, l->r), make_tuple(r->l, r->r)));
    }else if(l->left){
        __wspd(l->left, r);
        __wspd(l->right, r);
    }
}
