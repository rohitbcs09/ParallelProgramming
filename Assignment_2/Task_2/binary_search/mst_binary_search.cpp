#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <chrono>
#include <sstream>
#include <string>
#include <fstream>
#include <unordered_map>

using namespace std;

typedef struct edge {
    uint64_t u;
    uint64_t v;
    double w;
} Edge;

typedef std::vector<Edge*> EdgeList;

uint64_t g_seed = time(0);

uint64_t fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
}

Edge* createEdge(uint64_t u, uint64_t v, double w) {
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
    for(uint64_t i = 0; i<orig.size(); ++i){
        copy.push_back(copy_edge(orig[i]));
    }
}

void print(std::vector<Edge*> &edge_list) {
    for(uint64_t i = 0; i< edge_list.size(); ++i){
        std::cout << "( " << edge_list[i]->u << ", "
                  << edge_list[i]->v <<  ", " 
                  << edge_list[i]->w << " )";
    }
}

struct Comparator {
    bool operator()(Edge* &lhs, Edge* &rhs) {
        return lhs->w < rhs->w;
    }
};

void Par_Simulate_Priority_CW_BS(uint64_t n, std::vector<Edge*> &E, std::vector<uint64_t> &R) {
    std::vector<uint64_t> B(n+1, 0);
    std::vector<uint64_t> l(n+1, 0);
    std::vector<uint64_t> h(n+1, 0);

    std::vector<uint64_t> lo(n+1, 0);
    std::vector<uint64_t> hi(n+1, 0);
    std::vector<uint64_t> md(n+1, 0);

    uint64_t m = E.size();

    #pragma cilk grainsize = 8
    cilk_for(uint64_t u = 1; u<=n; ++u) {
        l[u] = 0;
        h[u] = m-1;
    }

    for(uint64_t k = 0; k< 1 + log2(m); ++k){
        #pragma cilk grainsize = 8
        cilk_for(uint64_t u = 1; u<=n; ++u) {
            B[u] = 0;
            lo[u] = l[u];
            hi[u] = h[u];
        }
        #pragma cilk grainsize = 8
        cilk_for(uint64_t i = 0; i< m; ++i) {
            uint64_t u = E[i]->u;
            md[u] = (lo[u] + hi[u])/2;
            if( i>=lo[u] && i<=md[u]) {
                B[u] = 1;
            }
        }
        #pragma cilk grainsize = 8
        cilk_for(uint64_t i = 0; i< m; ++i) {
            uint64_t u = E[i]->u;
            md[u] = (lo[u] + hi[u])/2;
            if( B[u] == 1 && i>=lo[u] && i<=md[u]) {
                h[u] = md[u];
            }
            else if(B[u] == 0 && i>=md[u] &&
                    i<=hi[u] ) {
                l[u] = md[u] + 1;
            }
        }
        #pragma cilk grainsize = 8
        cilk_for(uint64_t i = 0; i<m; ++i) {
            uint64_t u = E[i]->u;
            if(i==l[u]) {
                R[u] = i;
            }
        }
    }
}

void parallel_prefix_sum(std::vector<uint64_t> &arr, uint64_t nums, std::vector<uint64_t> &indexes) {
  
   if (nums == 1) {
       indexes[0] = arr[0];
       return;
   }

   std::vector<uint64_t> y(nums/2, 0);
   std::vector<uint64_t> z(nums/2, 0);

   cilk_for(int i = 0; i < nums/2; ++i) {
      y[i] = arr[2 * i] + arr[(2 * i) + 1];
   }

   parallel_prefix_sum(y, nums/2, z);

   cilk_for (int i = 0; i < nums; ++i) {
      if (i == 0) {
          indexes[0] = arr[0];
      } else if (i % 2 == 1) {
          indexes[i] = z[i / 2];
      } else {
          indexes[i] = z[(i - 1)/2] + arr[i];
      }
   }
}

void Par_Mst_Priority_CW(uint64_t n, std::vector<Edge*> &E, std::vector<uint64_t> &Mst) {
    std::vector<uint64_t> L(n+1, -1);
    std::vector<uint64_t> C(n+1, -1);
    std::vector<uint64_t> R(n+1, -1);

    #pragma cilk grainsize = 8
    cilk_for(uint64_t v = 1; v<=n; ++v) {
        L[v] = v;
    }
    uint64_t m = E.size();
    bool flag = m > 0 ? true : false;
    while(flag) {
        #pragma cilk grainsize = 8
        cilk_for(uint64_t v = 1; v<=n; ++v) {
           C[v] = fastrand() % 2;
        }

        Par_Simulate_Priority_CW_BS(n, E, R);
        #pragma cilk grainsize = 8
        cilk_for(uint64_t i = 0; i<m; ++i) {
           uint64_t u = E[i]->u; 
           uint64_t v = E[i]->v; 
           uint64_t w = E[i]->w; 
           if( C[u] == 0 && C[v] == 1 &&
               R[u] == i) {
               L[u] = v;
               Mst[i] = 1;
           }
        }
 
        #pragma cilk grainsize = 8
        cilk_for(uint64_t i = 0; i<m; ++i) {
            int u = E[i]->u;
            int v = E[i]->v;
            E[i]->u = L[u];
            E[i]->v = L[v];
            if(E[i]->u == E[i]->v) {
                E[i]->u = 0;
                E[i]->v = 0;
            }
        }
        flag = false;
        #pragma cilk grainsize = 8
        cilk_for(uint64_t i = 0; i< m; ++i) {
            if(E[i]->u != E[i]->v) {
                //std::cout << "IF: "<< i << " "<< E[i]->u << " " << E[i]->v << "\n";
                flag = true;
            }
        }
    }

    return;
}

int main(int argc, char** argv) {

    // Setting number of worker threads
    if(argc > 1) {
        __cilkrts_set_param("nworkers", argv[1]);
        std::cout << argv[1] << "\n";
    }

    EdgeList edge_list;

    std::string line;
    // std::ifstream infile("../input_graphs/as-skitter-in.txt");
    //std::ifstream infile("../input_graphs/com-amazon-in.txt");
    //std::ifstream infile("../input_graphs/com-friendster-in.txt");
    std::ifstream infile("../input_graphs/com-youtube-in.txt");
    //std::ifstream infile("temp_1.txt");
    std::getline(infile, line);
    std::istringstream iss(line);
    uint64_t n, m1;
    iss >> n >> m1;
    std::cout << n << " " << m1 << "\n";
    while(std::getline(infile, line)) {
        uint64_t u, v;
        double w;
        std::istringstream iss(line);
        iss >> u >> v >> w; 
        edge_list.push_back(createEdge(u, v, w));
        edge_list.push_back(createEdge(v, u, w));
    }
    uint64_t m = edge_list.size();
    std::vector<uint64_t> Mst(m, 0);
    
    std::sort(edge_list.begin(), edge_list.end(), Comparator());
    EdgeList copy_edge_list;
    deep_copy(copy_edge_list, edge_list);

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    Par_Mst_Priority_CW(n, edge_list, Mst);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "Running Time: " << time_span.count() << " seconds.\n";

    std::ofstream outfile;
    //outfile.open("output_mst_graphs/as-skitter-out.txt");
    //outfile.open("output_mst_graphs/com-amazon-in.txt");
    //outfile.open("output_mst_graphs/com-friendster-in.txt");
    outfile.open("output_mst_graphs/com-youtube-in.txt");
    //outfile.open("output_mst_graphs/temp_1.txt");
    outfile << "Running Time: " << time_span.count() << " seconds.\n";
    for(uint64_t i = 0; i<Mst.size(); ++i) {
        if(Mst[i]) {
            /*std::cout<< copy_edge_list[i]->u << " "
                     << copy_edge_list[i]->v << " "
                     << copy_edge_list[i]->w << "\n";*/
	    outfile << copy_edge_list[i]->u;
            outfile << " ";
	    outfile << copy_edge_list[i]->v;
            outfile << " ";
	    outfile << copy_edge_list[i]->w;
            outfile << " ";
	    //outfile << i;
	    outfile << "\n";
        }
    }
    outfile.close();
    return 1;
}
