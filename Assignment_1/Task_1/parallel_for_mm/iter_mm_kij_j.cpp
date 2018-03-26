#include<iostream>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#include<cilk/cilk.h>
#include<cilk/cilk_api.h>
#include <chrono>
#include <cstdlib>


int main(int argc, char *argv[]) {
	
	__cilkrts_set_param("nworkers", "68");
	int msize = atoi(argv[1]);
	
	int **a, **b, **res;
        a = new int * [msize];
        b = new int * [msize];
        res = new int * [msize];

        for (int ind = 0; ind < msize; ++ind) {
            a[ind] = new int [msize];
            b[ind] = new int [msize];
            res[ind] = new int [msize];
        }

	// make the two matrices.
	for (int i = 0; i < msize; i++ ) {
                for (int j = 0; j < msize; j++) {
                        a[i][j] = rand() % INT_MAX;
                        b[i][j] = rand() % INT_MAX;
			res[i][j] = 0;
                }
        }

	using namespace std::chrono;
    	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	//multiple the matrices
	for ( int k = 0; k < msize; k++) {
		for ( int i = 0; i < msize; i++) {
			cilk_for ( int j = 0; j < msize; j++) {
				res[i][j] = res[i][j] + a[i][k] * b[k][j];
			}
		}
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	
	std::cout << "CPU time used for kij : " << time_span.count() << std::endl;
	return 0;
}
