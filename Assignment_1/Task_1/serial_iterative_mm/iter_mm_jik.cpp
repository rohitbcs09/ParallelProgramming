#include<iostream>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#include <chrono>

#define msize 1024


int main() {

	int i, j, k;
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
	for (i = 0; i < msize; i++ ) {
                for (j = 0; j < msize; j++) {
                        a[i][j] = rand() % INT_MAX;
                        b[i][j] = rand() % INT_MAX;
			res[i][j] = 0;
                }
        }

	using namespace std::chrono;

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

	//multiple the matrices
	for ( j = 0; j < msize; j++) {
		for ( i = 0; i < msize; i++) {
			for ( k = 0; k < msize; k++) {
				res[i][j] = res[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	std::cout << "CPU time used for jik : " << time_span.count() << std::endl;
	return 0;
}