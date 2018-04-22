#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits.h>
#include <time.h>
#include <math.h>
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

int partition(std::vector<int> &arr, int left, int right) {
    
    int nums = right - left + 1;

    std::vector<int> res(nums, 0);
    
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
        if (arr[left + ind] < arr[right]) {
	    less_than_flag[ind] = 1;
	    equal_to_flag[ind] = 0;
	    greater_than_flag[ind] = 0;
	} else if (arr[left + ind] == arr[right]) {
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

void insertion_sort(std::vector<int> &arr, int start, int end) {
    int size = end - start + 1;
    if (size <= 1) {
        return;    
    }

    int i = start;
    while (i <= end) {
        int j = i;
        while (j >= 0 && arr[j - 1] > arr[j]) {
	    int temp = arr[j - 1];
            arr[j - 1] = arr[j];
            arr[j] = temp;
            j--;
        }
        i++;
    }

}


void quick_sort(std::vector<int> &arr, int left, int right, int serial_break) {
    if (left < right) {
        int nums = right - left + 1;
        if (nums > serial_break) {
            int rand_index = left + (fastrand() % nums);
            // Replace random index element with last element
            int temp = arr[rand_index];
            arr[rand_index] = arr[right];
            arr[right] = temp;
            int ind = partition(arr, left, right);
            cilk_spawn quick_sort(arr, left, ind - 1, serial_break);
            quick_sort(arr, ind + 1, right, serial_break);
        } else {
	    insertion_sort(arr, left, right);
        }
    }
}


int main(int argc, char** argv) {

    __cilkrts_set_param("nworkers", argv[1]);

    int nums = atoi(argv[2]);
    int serial_break = atoi(argv[3]);

    std::vector<int> input(nums, 0);
    fill_input(input, nums);
    

    //print_arr(input, nums);
    
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now(); 

    quick_sort(input, 0, nums - 1, serial_break);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

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
