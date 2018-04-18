#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cilk/cilk.h>
#include <chrono>

using namespace std;

typedef struct edge {
    int u;
    int v;
    int w;
} Edge;

typedef std::vector<Edge*> EdgeList;

uint64_t g_seed = time(0);

uint64_t fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
}

Edge* createEdge(int u, int v, int w) {
    Edge *edge = new Edge;
    edge->u = u;
    edge->v = v;
    edge->w = w;
    return edge;
}

Edge* copy_edge(Edge *orig) {
    Edge *copy = new Edge;
    copy->u = orig->u;
    copy->v = orig->v;
    copy->w = orig->w;
    return copy;
}

void deep_copy(std::vector<Edge*> &copy, 
               std::vector<Edge*> &orig) {          
    for(int i = 0; i<orig.size(); ++i){
        copy.push_back(copy_edge(orig[i]));
    }
}

void print(std::vector<Edge*> &edge_list) {
    for(int i = 0; i< edge_list.size(); ++i){
        std::cout << "( " << edge_list[i]->u << ", "
                  << edge_list[i]->v <<  ", " 
                  << edge_list[i]->w << " )\n";
    }
}

struct Comparator {
    bool operator()(Edge* &lhs, Edge* &rhs) {
        return lhs->w < rhs->w;
    }
};

void Par_Simulate_Priority_CW_BS(int n, std::vector<Edge*> &E, std::vector<int> &R) {
    std::vector<int> B(n, 0);
    std::vector<int> l(n, 0);
    std::vector<int> h(n, 0);

    std::vector<int> lo(n, 0);
    std::vector<int> hi(n, 0);
    std::vector<int> md(n, 0);

    int m = E.size();

    cilk_for(int u = 0; u<n; ++u) {
        l[u] = 1;
        h[u] = m;
    }

    for(int k = 1; k< 1 + log2(m); ++k){
        cilk_for(int u = 0; u<n; ++u) {
            B[u] = 0;
            lo[u] = l[u];
            hi[u] = h[u];
        }
        cilk_for(int i = 0; i< m; ++i) {
            int u = E[i]->u;
            md[u] = (lo[u] + hi[u])/2;
            if( i>=lo[u] && i<=md[u]) {
                B[u] = 1;
            }
        }
        cilk_for(int i = 0; i< m; ++i) {
            int u = E[i]->u;
            md[u] = (lo[u] + hi[u])/2;
            if( B[u] == 1 && i>=lo[u] && i<=md[u]) {
                h[u] = md[u];
            }
            else if(B[u] == 0 && i>=md[u] &&
                    i<=hi[u] ) {
                l[u] = md[u] + 1;
            }
        }
        cilk_for(int i = 0; i<m; ++i) {
            int u = E[i]->u;
            if(i==l[u]) {
                R[u] = i;
            }
        }
    }
}


void Par_Mst_Priority_CW(int n, std::vector<Edge*> E, std::vector<int> &Mst) {
    std::vector<int> L(n, -1);
    std::vector<int> C(n, -1);
    std::vector<int> R(n, -1);

    //std::sort(E.begin(), E.end(), Comparator());
    cilk_for(int v = 0; v<n; ++v) {
        L[v] = v;
    }
    int m = E.size();
    bool flag = m > 0 ? true : false;
    while(flag) {
        cilk_for(int v = 0; v<n; ++v) {
           C[v] = fastrand() % 2;
        }

        Par_Simulate_Priority_CW_BS(n, E, R);

        cilk_for(int i = 0; i<m; ++i) {
           int u = E[i]->u; 
           int v = E[i]->v; 
           int w = E[i]->w; 
           if( C[u] == 0 && C[v] == 1 &&
               R[u] == i) {
               L[u] = v;
               Mst[i] = 1;
           }
        }
        cilk_for(int i = 0; i<m; ++i) {
            E[i]->u = L[E[i]->u];
            E[i]->v = L[E[i]->v];
        }
        flag = false;
        cilk_for(int i = 0; i< m; ++i) {
            if(E[i]->u != E[i]->v) {
                flag = true;
            }
        }
    }
    return;
}

int main() {

    EdgeList edge_list;

    edge_list.push_back(createEdge(0, 1, 2));
    edge_list.push_back(createEdge(1, 0, 2));

    edge_list.push_back(createEdge(1, 2, 3));
    edge_list.push_back(createEdge(2, 1, 3));

    edge_list.push_back(createEdge(0, 3, 6));
    edge_list.push_back(createEdge(3, 0, 6));

    edge_list.push_back(createEdge(1, 3, 8));
    edge_list.push_back(createEdge(3, 1, 8));

    edge_list.push_back(createEdge(3, 4, 9));
    edge_list.push_back(createEdge(4, 3, 9));

    edge_list.push_back(createEdge(1, 4, 5));
    edge_list.push_back(createEdge(4, 1, 5));

    edge_list.push_back(createEdge(2, 4, 7));
    edge_list.push_back(createEdge(4, 2, 7));
    
    int m = edge_list.size();
    std::vector<int> Mst(m, 0);
    
    //print(edge_list);
    std::sort(edge_list.begin(), edge_list.end(), Comparator());
    EdgeList copy_edge_list;
    deep_copy(copy_edge_list, edge_list);

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    Par_Mst_Priority_CW(5, edge_list, Mst);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "Running Time: " << time_span.count() << " seconds.\n";

    for(int i = 0; i<Mst.size(); ++i) {
        if(Mst[i]) {
            std::cout<< copy_edge_list[i]->u << " "
                     << copy_edge_list[i]->v << " "
                     << copy_edge_list[i]->w << "\n";
        }
    }
    return 1;
}
