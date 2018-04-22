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
#include <cmath>

using namespace std;

uint64_t toss_count = 0;

uint64_t processor = 1;

typedef struct edge {
    uint64_t u;
    uint64_t v;
    double w;
} Edge;

typedef std::vector<Edge*> EdgeList;

uint64_t g_seed;

uint64_t fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
}

void Par_Counting_Rank(std::vector<uint64_t> &S, uint64_t nums, uint64_t d,
                        std::vector<uint64_t> &r, uint64_t processor);


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

uint64_t EXTRACT_BIT_SEGMENT(uint64_t num, uint64_t start_bit, uint64_t end_bit) {
    uint64_t mask = ~(~0 << (end_bit - start_bit + 1));
    return mask & (num >> start_bit);
}

void Par_radix_sort(std::vector<uint64_t> &arr, uint64_t nums , uint64_t bits, uint64_t processor) {
    std::vector<uint64_t> S(nums);
    std::vector<uint64_t> r(nums);
    std::vector<uint64_t> B(nums);

    uint64_t d = ceil(log((nums * 1.0) / processor * log (nums)));

    for(uint64_t k = 0; k < bits; ++k) {
        uint64_t q = (k + d <= bits) ? d : (bits - k);

        //#pragma cilk grainsize = 2048
        for(uint64_t i = 0; i < nums; ++i) {
	    S[i] = EXTRACT_BIT_SEGMENT(arr[i], k, k + q - 1);    
        }

 	Par_Counting_Rank(S, nums, q, r, processor );

        //#pragma cilk grainsize = 2048
        for (uint64_t i = 0; i < nums; ++i) {
	    B[r[i]] = arr[i];
        }

        //#pragma cilk grainsize = 2048
        for(uint64_t i = 0; i < nums; ++i) {
	    arr[i] = B[i];	
        }
    }

}

void Par_Simulate_Priority_CW_RS_2(uint64_t n, std::vector<Edge*> &E,
         std::vector<uint64_t> &R, uint64_t processor) {

    uint64_t m = E.size();
    std::vector<uint64_t> A(m);
    uint64_t k = std::ceil(log2(m)) + 1;

    #pragma cilk grainsize = 2048
    cilk_for(uint64_t i = 0; i<m; ++i) {
        A[i] = (E[i]->u << k) + i;
    }
    Par_radix_sort(A, m, k + std::ceil(log2(n)), processor); 
    
    #pragma cilk grainsize = 2048
    cilk_for(uint64_t i = 0; i<m; ++i) {
        uint64_t u = (A[i] >> k);
        uint64_t j = A[i] - (u<<k);
        if (i == 0 || u != (A[i-1] >> k)) {
            R[u] = j;
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

   //#pragma cilk grainsize = 2048
   for(uint64_t i = 0; i < nums/2; ++i) {
      y[i] = arr[2 * i] + arr[(2 * i) + 1];
   }

   parallel_prefix_sum(y, nums/2, z);

   //#pragma cilk grainsize = 2048
   for (uint64_t i = 0; i < nums; ++i) {
      if (i == 0) {
          indexes[0] = arr[0];
      } else if (i % 2 == 1) {
          indexes[i] = z[i / 2];
      } else {
          indexes[i] = z[(i - 1)/2] + arr[i];
      }
   }
}


void Par_Counting_Rank(std::vector<uint64_t> &S, uint64_t nums, uint64_t d,
                        std::vector<uint64_t> &r, uint64_t processor) {

    std::vector<std::vector<uint64_t>> f(pow(2, d), 
                std::vector<uint64_t> (processor));
    std::vector<std::vector<uint64_t>> r_1(pow(2, d),
                std::vector<uint64_t> (processor));   
    std::vector<uint64_t> js(processor);
    std::vector<uint64_t> je(processor);
    std::vector<uint64_t> ofs(processor);

    //#pragma cilk grainsize = 2048
    for(uint64_t i=0; i<processor; ++i) {
        for (uint64_t j = 0; j<pow(2,d); ++j) {
            f[j][i] = 0;
        }
        js[i] = i * (floor(nums / processor));
        je[i] = i < processor ?  (i+1) * (floor(nums / processor)) - 1 : nums - 1;

        for (uint64_t j = js[i]; j <= je[i]; ++j) {
            f[S[j]][i] = f[S[j]][i] + 1;
        }
    }

    for (uint64_t j = 0; j < pow(2, d); ++j) {
            std::vector<uint64_t> temp(processor);
            parallel_prefix_sum(f[j], processor, temp);
            f[j] = temp;
    }
     
    
    //#pragma cilk grainsize = 2048
    for (uint64_t i = 0; i < processor; ++i) {
        ofs[i] = 0;
        for (uint64_t j = 0; j < pow(2, d); ++j) {
           r_1[j][i] = (i == 0) ? ofs[i] : ofs[i] + f[j][i - 1];
            ofs[i] = ofs[i] + f[j][processor - 1];
        }
        
        for (uint64_t j = js[i]; j <= je[i]; ++j) {
            r[j] = r_1[S[j]][i];
            r_1[S[j]][i] = r_1[S[j]][i] + 1;
        }
    }
}


void Par_Mst_Priority_CW(uint64_t n, std::vector<Edge*> &E,
         std::vector<uint64_t> &Mst, uint64_t processor) {

    std::vector<uint64_t> L(n+1);
    std::vector<uint64_t> C(n+1);
    std::vector<uint64_t> R(n+1);

    #pragma cilk grainsize = 2048
    cilk_for(uint64_t v = 1; v<=n; ++v) {
        L[v] = v;
    }
    uint64_t m = E.size();
    bool flag = m > 0 ? true : false;
    while(flag) {
        #pragma cilk grainsize = 2048
        cilk_for(uint64_t v = 1; v<=n; ++v) {
           C[v] = fastrand() % 2;
        }

        Par_Simulate_Priority_CW_RS_2(n, E, R, processor);
        #pragma cilk grainsize = 2048
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
 
        #pragma cilk grainsize = 2048
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
        #pragma cilk grainsize = 2048
        cilk_for(uint64_t i = 0; i< m; ++i) {
            if(E[i]->u != E[i]->v) {
                flag = true;
            }
        }
        ++toss_count;
    }
    return;
}

int main(int argc, char** argv) {
    
    // Setting number of worker threads
    if(argc > 1) {
        __cilkrts_set_param("nworkers", argv[1]);
        processor = atoi(argv[1]);        
        std::cout << "Num Processors: "  << processor << " \n";
    }

    EdgeList edge_list;

    std::string line;
    std::ifstream infile("../input_graphs/com-youtube-in.txt");
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
    srand(time(NULL));
    g_seed=rand();
    uint64_t m = edge_list.size();
    std::vector<uint64_t> Mst(m, 0);
    
    std::sort(edge_list.begin(), edge_list.end(), Comparator());
    EdgeList copy_edge_list;
    deep_copy(copy_edge_list, edge_list);

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    Par_Mst_Priority_CW(n, edge_list, Mst, processor);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    std::cout << "Running Time: " << time_span.count() << " seconds.\n";

    std::ofstream outfile;
    outfile.open("output_mst_graphs/com-youtube-in.txt");
    outfile << "Running Time: " << time_span.count() << " seconds.\n";
    for(uint64_t i = 0; i<Mst.size(); ++i) {
        if(Mst[i]) {
            outfile << copy_edge_list[i]->u;
            outfile << " ";
            outfile << copy_edge_list[i]->v;
            outfile << " ";
            outfile << copy_edge_list[i]->w;
            outfile << " ";
            outfile << "\n";
        }
    }
    std::cout << "Toss Count = " << toss_count ;
    outfile.close();
    return 1;
}