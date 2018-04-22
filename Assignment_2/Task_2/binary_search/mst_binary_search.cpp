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

typedef struct edge {
    uint64_t u;
    uint64_t v;
    double w;
} Edge;

typedef std::vector<Edge*> EdgeList;

uint64_t g_seed;

int g_index;

uint64_t fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
}

int fastrand_index() {
  g_index = (214013 * g_index + 2531011); 
  return (g_index>>16) & 0x7FFF; 
}

void fill_input(std::vector<int> &arr, int msize) {
   for(int i = 0; i < msize; ++i) {
       arr[i] = fastrand_index();
   }
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

void Par_Simulate_Priority_CW_BS(uint64_t n, std::vector<Edge*> &E, std::vector<uint64_t> &R) {
    std::vector<uint64_t> B(n+1, 0);
    std::vector<uint64_t> l(n+1, 0);
    std::vector<uint64_t> h(n+1, 0);

    std::vector<uint64_t> lo(n+1, 0);
    std::vector<uint64_t> hi(n+1, 0);
    std::vector<uint64_t> md(n+1, 0);

    uint64_t m = E.size();

    #pragma cilk grainsize = 2048
    cilk_for(uint64_t u = 1; u<=n; ++u) {
        l[u] = 0;
        h[u] = m-1;
    }

    for(uint64_t k = 0; k< 1 + log2(m); ++k){
        #pragma cilk grainsize = 2048
        cilk_for(uint64_t u = 1; u<=n; ++u) {
            B[u] = 0;
            lo[u] = l[u];
            hi[u] = h[u];
        }
        #pragma cilk grainsize = 2048
        cilk_for(uint64_t i = 0; i< m; ++i) {
            uint64_t u = E[i]->u;
            md[u] = (lo[u] + hi[u])/2;
            if( i>=lo[u] && i<=md[u]) {
                B[u] = 1;
            }
        }
        #pragma cilk grainsize = 2048
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
        #pragma cilk grainsize = 2048
        cilk_for(uint64_t i = 0; i<m; ++i) {
            uint64_t u = E[i]->u;
            if(i==l[u]) {
                R[u] = i;
            }
        }
    }
}

void parallel_prefix_sum(std::vector<int> &arr, int nums, std::vector<int> &indexes) {
  
   /* 
   int num = 0;
   for (int ind = 0; ind < nums; ++ind) {
       num = num + arr[ind];
       indexes[ind] = num;
   } */

   if (nums == 1) {
       indexes[0] = arr[0];
       return;
   }

   std::vector<int> y(nums/2, 0);
   std::vector<int> z(nums/2, 0);

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

int partition(std::vector<Edge*> &arr, int left, int right) {
    
    int nums = right - left + 1;

    std::vector<Edge*> res(nums);
    
    std::vector<int> less_than_flag(nums, 0);
    std::vector<int> equal_to_flag(nums, 0);
    std::vector<int> greater_than_flag(nums, 0);

    std::vector<int> less_than_index(nums, 0);
    std::vector<int> equal_to_index(nums, 0);
    std::vector<int> greater_than_index(nums, 0);

    cilk_for (int ind = 0; ind < nums; ind++) {
        res[ind] = arr[left + ind];
    }

    cilk_for (int ind = 0; ind < nums; ++ind) {
        if (arr[left + ind]->w < arr[right]->w) {
	    less_than_flag[ind] = 1;
	    equal_to_flag[ind] = 0;
	    greater_than_flag[ind] = 0;
	} else if (arr[left + ind]->w == arr[right]->w) {
	    less_than_flag[ind] = 0;
            equal_to_flag[ind] = 1;
            greater_than_flag[ind] = 0;
        } else {
	    less_than_flag[ind] = 0;
            equal_to_flag[ind] = 0;
            greater_than_flag[ind] = 1;
        }
    }
   parallel_prefix_sum(less_than_flag, nums, less_than_index);
   parallel_prefix_sum(equal_to_flag, nums, equal_to_index);
   parallel_prefix_sum(greater_than_flag, nums, greater_than_index);

   int lt_index_max = less_than_index[nums - 1];
   arr[left + lt_index_max] = res[nums - 1];

    cilk_for (int ind = 0; ind < nums; ind++) {
        if (less_than_flag[ind]) {
            arr[left + less_than_index[ind] - 1] = res[ind];
        } else if (equal_to_flag[ind]) {
            arr[left + lt_index_max + equal_to_index[ind] - 1] = res[ind];
        } else if (greater_than_index[ind]) {
            arr[left + greater_than_index[ind] + lt_index_max + equal_to_index[nums - 1] - 1] = res[ind];
        }
    }

    return left + lt_index_max;
}

void insertion_sort(std::vector<Edge*> &arr, int start, int end) {
    int size = end - start + 1;
    if (size <= 1) {
        return;    
    }
    int i = start;
    while (i <= end) {
        int j = i;
        while (j > start && arr[j - 1]->w > arr[j]->w) {
	    Edge* temp = arr[j - 1];
            arr[j - 1] = arr[j];
            arr[j] = temp;
            j--;
        }
        i++;
    }
}

void quick_sort(std::vector<Edge*> &arr, int left, int right, int serial_break) {
    if (left < right) {
        int nums = right - left + 1;
        if (nums > serial_break) {
            //int rand_index = left + (fastrand_index() % nums);
            // Replace random index element with last element
            //Edge *temp = arr[rand_index];
            //arr[rand_index] = arr[right];
            //arr[right] = temp;
            int ind = partition(arr, left, right);
            /*cilk_spawn*/ quick_sort(arr, left, ind - 1, serial_break);
            quick_sort(arr, ind + 1, right, serial_break);
        } else {
	    insertion_sort(arr, left, right);
        }
    }
}

void Par_Mst_Priority_CW(uint64_t n, std::vector<Edge*> &E, std::vector<uint64_t> &Mst) {
    std::vector<uint64_t> L(n+1, -1);
    std::vector<uint64_t> C(n+1, -1);
    std::vector<uint64_t> R(n+1, -1);

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

        Par_Simulate_Priority_CW_BS(n, E, R);
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
    }

    return;
}

int main(int argc, char** argv) {

    // Setting number of worker threads
    if(argc > 1) {
        __cilkrts_set_param("nworkers", argv[1]);
        std::cout << argv[1] << "\n";
    }

    if(argc < 4) {
        std::cout << "Missing Input params - <binary> <processor> <input_file> <out_file>\n";
        return 1;
    }

    EdgeList edge_list;
    EdgeList copy_edge_list;

    std::string line;
    std::ifstream infile(argv[2]);
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
        copy_edge_list.push_back(createEdge(u, v, w));
        copy_edge_list.push_back(createEdge(v, u, w));
    }
    uint64_t m = edge_list.size();
    std::vector<uint64_t> Mst(m, 0);

    srand(time(NULL));
    g_seed=rand();
    g_index=time(0);

    quick_sort(edge_list, 0, edge_list.size() - 1, 1024);

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    Par_Mst_Priority_CW(n, edge_list, Mst);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "Running Time: " << time_span.count() << " seconds.\n";

    std::ofstream outfile;
    outfile.open(argv[3]);
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
    outfile.close();
    return 1;
}
