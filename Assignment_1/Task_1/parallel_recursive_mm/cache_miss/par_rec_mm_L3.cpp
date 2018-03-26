#include <iostream>
#include <vector>
#include <ctime>
#include <stdint.h>
#include <time.h>
#include <cstdlib>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <papi.h>
#include <chrono>

using namespace std;

void handle_error(int err){
    std::cerr << "PAPI error: " << err << std::endl;
}

/*
FAST RANDOM NUMBER GENERATOR:
    o https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
*/


uint64_t g_seed = time(0);

static inline uint64_t fastrand() { 
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
} 

void fillMatrix(std::vector<std::vector<uint64_t> > &arr, int msize) {
    for(int i = 0; i < msize; ++i) {
        for(int j = 0; j < msize; ++j) {
            arr[i][j] = j; //fastrand() % 10;
        }
    }
}

void print_Matrix(std::vector<std::vector<uint64_t> > &arr, int x_s, int x_e, 
                 int y_s, int y_e) {
     for(int i = x_s; i<x_e; ++i) {
        for(int j = y_s; j<y_e; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }   
}

void printMatrix(std::vector<std::vector<uint64_t> > &arr, int n) {
     std::cout << "\n";
     for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }   
}

void SUM_SUB_MATRIX(std::vector<std::vector<uint64_t> > &X, 
                    std::vector<std::vector<uint64_t> > &X_, int x_s, int x_e, 
                    int y_s, int y_e) {
    //print_Matrix(X, x_s, x_e, y_s, y_e);
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            X[i][j] += X_[row][col];
        }
    }
}

void COPY_SUB_MATRIX(std::vector<std::vector<uint64_t> > &X_,
                                   std::vector<std::vector<uint64_t> > &X, int x_s,
                                   int x_e, int y_s, int y_e) {
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            X_[row][col] = X[i][j];
        }
    }
}

static inline void SUM_MATRIX(std::vector<std::vector<uint64_t> > &Z,
                              std::vector<std::vector<uint64_t> > &T) {
    for(int i = 0; i< Z.size(); ++i) {
        for(int j = 0; j < Z[0].size(); ++j) {
            Z[i][j] = Z[i][j] + T[i][j];
        }
    }
}

void Matrix_Multiply(std::vector<std::vector<uint64_t> > &X, 
                     std::vector<std::vector<uint64_t> > &Y,
                     std::vector<std::vector<uint64_t> > &Z, int x_row, 
					 int x_col, int y_row, int y_col, int z_row, int z_col, 
					 int n) {

    for(int i = 0; i<n; ++i){
        for(int k = 0; k<n; ++k) {
            for(int j = 0; j<n; ++j) {
                Z[z_row + i][z_col + j] += X[x_row + i][x_col + k] * Y[y_row + k][y_col + j];
            }
        }
    }
}

void PAR_REC_MEM(std::vector<std::vector<uint64_t> > &X,
                 std::vector<std::vector<uint64_t> > &Y,
				 std::vector<std::vector<uint64_t> > &Z, int x_row, int x_col,
				 int y_row, int y_col, int z_row, int z_col, int n, int base) {
        
    if(n == base) {
        Matrix_Multiply(X, Y, Z, x_row, x_col, y_row, y_col, z_row, z_col, base);
    }
    else {
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row, x_col, y_row, y_col, z_row, z_col, n/2, base);
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row, x_col, y_row, y_col + n/2, z_row, z_col + n/2, n/2, base);
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col, y_row, y_col, z_row + n/2, z_col, n/2, base);
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col, y_row, y_col + n/2, z_row + n/2, z_col + n/2, n/2, base);
        cilk_sync;
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row, x_col + n/2, y_row + n/2, y_col, z_row, z_col, n/2, base);
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row, x_col + n/2, y_row + n/2, y_col + n/2, z_row, z_col + n/2, n/2, base);
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col + n/2, y_row + n/2, y_col, z_row + n/2, z_col, n/2, base);
        cilk_spawn PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col + n/2, y_row + n/2, y_col + n/2, z_row + n/2, z_col + n/2, n/2, base);
        cilk_sync;
    }
}

int main(int argc, char *argv[]) {

    __cilkrts_set_param("nworkers", "68");

    int msize = atoi(argv[1]);
    int base = atoi(argv[2]);

    std::vector< std::vector<uint64_t> > A(msize, std::vector<uint64_t>(msize, 0));
    std::vector< std::vector<uint64_t> > B(msize, std::vector<uint64_t>(msize, 0));
    std::vector< std::vector<uint64_t> > Z(msize, std::vector<uint64_t>(msize, 0));
    
    fillMatrix(A, msize);
    fillMatrix(B, msize);

    // Print array
    // printMatrix(A, 8);
    // printMatrix(B, 8);

    int numEvents = 1;
    long long values[1];
    int events[1] = {PAPI_L3_DCM};

    using namespace std::chrono;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
 
    if (PAPI_start_counters(events, numEvents) != PAPI_OK) {
	 // std::cout << PAPI_start_counters(events, numEvents) << std::endl;
         handle_error(1);
     }
 
    
    cilk_spawn PAR_REC_MEM(A, B, Z, 0, 0, 0, 0, 0, 0, msize, base);
    cilk_sync;
    // printMatrix(Z, 8);

    if ( PAPI_stop_counters(values, numEvents) != PAPI_OK) {
            handle_error(1);
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    std::cout << "CPU time used : " << time_span.count() << std::endl;
    std::cout<<"L3 misses: "<< values[0] << std::endl;

    return 1;
}
