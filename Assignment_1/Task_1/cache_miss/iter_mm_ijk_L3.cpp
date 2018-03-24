#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <papi.h>

#define msize 1024

void handle_error(int err){
    std::cerr << "PAPI error: " << err << std::endl;
}


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

	int numEvents = 1;
    	long long values[1];
    	int events[1] = {PAPI_L3_TCM};

	if (PAPI_start_counters(events, numEvents) != PAPI_OK) {
            handle_error(1);
	}

	//multiple the matrices
	start = clock();
	for ( i = 0; i < msize; i++) {
		for ( j = 0; j < msize; j++) {
			for ( k = 0; k < msize; k++) {
				res[i][j] = res[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	end = clock();

	if ( PAPI_stop_counters(values, numEvents) != PAPI_OK) {
            handle_error(1);
	}

        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used for ijk : " << cpu_time_used << std::endl;
	std::cout<<"L3 misses: "<<values[0]<<std::endl;

	return 0;
}
