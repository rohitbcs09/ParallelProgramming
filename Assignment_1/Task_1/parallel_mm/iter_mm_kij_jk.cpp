#include<iostream>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#include<cilk/cilk.h>

#define msize 1024


int main() {

	int **a, **b, **res;
        a = new int * [msize];
        b = new int * [msize];
        res = new int * [msize];

        for (int ind = 0; ind < msize; ++ind) {
            a[ind] = new int [msize];
            b[ind] = new int [msize];
            res[ind] = new int [msize];
        }

	clock_t start, end;
 	double cpu_time_used;
	// make the two matrices.
	for (int i = 0; i < msize; i++ ) {
                for (int j = 0; j < msize; j++) {
                        a[i][j] = rand() % INT_MAX;
                        b[i][j] = rand() % INT_MAX;
			res[i][j] = 0;
                }
        }

	//multiple the matrices
	start = clock();
	cilk_for ( int k = 0; k < msize; k++) {
		for ( int i = 0; i < msize; i++) {
			cilk_for ( int j = 0; j < msize; j++) {
				res[i][j] = res[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	end = clock();
	
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used for kij : " << cpu_time_used << std::endl;
	return 0;
}
