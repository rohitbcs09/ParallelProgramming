#include<iostream>
#include<stdlib.h>
#include<limits.h>
#include<time.h>

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

	clock_t start, end;
 	double cpu_time_used;
	// make the two matrices.
	for (i = 0; i < msize; i++ ) {
                for (j = 0; j < msize; j++) {
                        a[i][j] = rand() % INT_MAX;
                        b[i][j] = rand() % INT_MAX;
			res[i][j] = 0;
                }
        }


	//multiple the matrices
	start = clock();
	for ( j = 0; j < msize; j++) {
		for ( i = 0; i < msize; i++) {
			for ( k = 0; k < msize; k++) {
				res[i][j] = res[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	end = clock();
	

        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used for jik : " << cpu_time_used << std::endl;
	return 0;
}
