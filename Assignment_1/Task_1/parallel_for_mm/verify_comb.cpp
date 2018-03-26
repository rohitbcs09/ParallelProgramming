#include<iostream>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#include<cilk/cilk.h>

#define msize 512


int main() {

	int **a, **b, **res, **ver_res;
        a = new int * [msize];
        b = new int * [msize];
        res = new int * [msize];
	ver_res = new int * [msize];

        for (int ind = 0; ind < msize; ++ind) {
            a[ind] = new int [msize];
            b[ind] = new int [msize];
            res[ind] = new int [msize];
	    ver_res[ind] = new int [msize];
        }
	
	// make the two matrices.
	for (int i = 0; i < msize; i++ ) {
                for (int j = 0; j < msize; j++) {
                        a[i][j] = rand() % INT_MAX;
                        b[i][j] = rand() % INT_MAX;
			res[i][j] = 0;
			ver_res[i][j] = 0;
                }
        }


	//multiple the matrices
	cilk_for ( int i = 0; i < msize; i++) {
		for ( int k = 0; k < msize; k++) {
			for ( int j = 0; j < msize; j++) {
				res[i][j] = res[i][j] + a[i][k] * b[k][j];
			}
		}
	}


	for ( int i = 0; i < msize; i++) {
                for ( int k = 0; k < msize; k++) {
                        for ( int j = 0; j < msize; j++) {
                                ver_res[i][j] = ver_res[i][j] + a[i][k] * b[k][j];
                        }
                }
        }

	for (int i = 0; i < msize; i++) {
		for (int j = 0; j < msize; j++) {
			if ( res[i][j] != ver_res[i][j]) {
				std::cout << " NOT MAtching" << std::endl;
				return 0;
			}
		}
	}


	std::cout << " Matching"<< std::endl;
	return 0;
}
