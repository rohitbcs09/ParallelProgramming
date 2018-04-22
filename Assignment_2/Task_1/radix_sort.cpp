#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include<cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <chrono>

using namespace std;

uint64_t g_seed = time(0);

static inline uint64_t fastrand() { 
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
} 

void fill_input(std::vector<int> &arr, int msize) {
   for(int i = 0; i < msize; ++i) {
       arr[i] = fastrand() % 1024;
   }
}

void print_arr(std::vector<int> &arr, int msize) {
    std::cout << std::endl;
    for(int i = 0; i < msize; ++i) {
       std::cout << arr[i] << " ";
    }
    std::cout << "\n";
}

void parallel_prefix_sum(std::vector<int> &arr, int nums, std::vector<int> &indexes) {
  
   /* 
   int num = 0;
   for (int ind = 0; ind < nums; ++ind) {
       num = num + arr[ind];
       indexes[ind] = num;
   }*/ 

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

void Par_Counting_Rank ( std::vector<int> &S, int nums, int d, std::vector<int> &r, int processor ) {
    std::vector<std::vector<int>> f((int) pow(2, d), std::vector<int> (processor, 0));
    std::vector<std::vector<int>> r_1((int) pow(2, d), std::vector<int> (processor, 0));   

    std::vector<int> js(processor, 0);
    std::vector<int> je(processor, 0);
    std::vector<int> ofs(processor, 0);

    cilk_for (int i = 0; i < processor; ++i) {
        for (int j = 0; j < (int) pow(2, d); ++j) {
            f[j][i] = 0;
        }
        
	js[i] = i * ((int) floor(nums / processor));
        je[i] = i < processor - 1 ?  (i + 1) * ((int) floor(nums / processor)) - 1 : nums - 1;

        //std::cout << "start - " << js[i] << " end - " << je[i] << std::endl;
   

        for (int j = js[i]; j <= je[i]; ++j) {
	    f[S[j]][i] = f[S[j]][i] + 1;
        }
    }

    for (int j = 0; j < (int) pow(2, d); ++j) {
            std::vector<int> temp(processor, 0);
            parallel_prefix_sum(f[j], processor, temp);
	    //print_arr(temp, 3);
            f[j] = temp;
    }
     
    
    cilk_for (int i = 0; i < processor; ++i) {
        ofs[i] = 0;
        for (int j = 0; j < (int) pow(2, d); ++j) {
   	    r_1[j][i] = (i == 0) ? ofs[i] : ofs[i] + f[j][i - 1];
            ofs[i] = ofs[i] + f[j][processor - 1];
        }
        
        for (int j = js[i]; j <= je[i]; ++j) {
            r[j] = r_1[S[j]][i];
            r_1[S[j]][i] = r_1[S[j]][i] + 1;
        }
    }

}


int EXTRACT_BIT_SEGMENT(int num, int start_bit, int end_bit) {
    unsigned long mask = ~(~0 << (end_bit - start_bit + 1));
    int val = mask & (num >> start_bit);
    //std::cout << "num - " << num << " start_bit - " << start_bit << " - end_bit - " << end_bit << " - segment - " << val << std::endl;
    return mask & (num >> start_bit);
}

void radix_sort(std::vector<int> &arr, int nums , int bits, int processor) {
    std::vector<int> S(nums, 0);
    std::vector<int> r(nums, 0);
    std::vector<int> B(nums, 0);

    int d = ceil(log((nums * 1.0) / processor * log (nums)));

    for (int k = 0; k < bits; ++k) {
        int q = (k + d <= bits) ? d : (bits - k);

        cilk_for (int i = 0; i < nums; ++i) {
	    S[i] = EXTRACT_BIT_SEGMENT(arr[i], k, k + q - 1);    
        }

        //std::cout << k << " - " << k + q - 1<< std::endl;
        //print_arr(S, nums); 

 	Par_Counting_Rank ( S, nums, q , r, processor );
  
        //std::cout << "********rank***********" << std::endl;
        //print_arr(r, nums);
        //std::cout << "************************" << std::endl;
        
	cilk_for (int i = 0; i < nums; ++i) {
	    B[ r[ i ] ] = arr[ i ];
        }

        cilk_for (int i = 0; i < nums; ++i) {
	    arr[ i ] = B[ i ];	
        }
        //std::cout << "******* nums **************" << std::endl;
        //print_arr(arr, nums);
        //std::cout << "************************" << std::endl;
	
    }

}


int main(int argc, char** argv) {

    __cilkrts_set_param("nworkers", argv[1]);

    int nums = atoi(argv[2]);
    std::vector<int> input(nums, 0);
    fill_input(input, nums);
    
    //print_arr(input, nums);
    
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now(); 

    radix_sort(input, nums, 10, atoi(argv[1]));
    
    //std::vector<int> res(20, 0);
    //Par_Counting_Rank(input, 20, 10, res, atoi(argv[1]));

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    
    //print_arr(res, 20);
    //print_arr(input, nums);
    
    bool sorted = true;
    for (int i = 1; i < nums; ++i) {
      if (input[i - 1] > input[i]) {
          sorted = false;
          break;
      }
    } 

    if (sorted) {
        std::cout << "Array sorted " << std::endl;
    } else {
        std::cout << "Array unsorted " << std::endl;
    }
    //print_arr(input, nums);

    std::cout << "CPU time used for rand-quicksort : " << time_span.count() <<  " with core " << argv[1] << std::endl;
    return 0;
}
