#include "../include/rec_mat_mul.h"
#include <ctime>
#include <stdint.h>
using namespace std;

/*

FAST RANDOM NUMBER GENERATOR:
    o https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c

*/

#define ROWS 8 
#define COLS 8 

// Input Matrices
Matrix X(ROWS, std::vector<uint64_t>(COLS, 0));
Matrix Y(ROWS, std::vector<uint64_t>(COLS, 0));
Matrix Z(ROWS, std::vector<uint64_t>(COLS, 0));

// Thread Pool 
std::vector<Task *> pool;

uint64_t g_seed = time(0);

void PAR_REC_MEM(Matrix &X, Matrix &Y, Matrix &Z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, int id);
 
static inline uint64_t fastrand() { 
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
} 

void fillMatrix(Matrix &arr) {
    for(int i = 0; i<ROWS; ++i) {
        for(int j = 0; j<COLS; ++j) {
            arr[i][j] = j; //fastrand() % 10;
        }
    }
}

void print_Matrix(Matrix &arr, int x_s, int x_e, 
                 int y_s, int y_e) {
     for(int i = x_s; i<x_e; ++i) {
        for(int j = y_s; j<y_e; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }   
}

void printMatrix(Matrix &arr, int n) {
     std::cout << "\n";
     for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }   
}

void SUM_SUB_MATRIX(Matrix &X, 
                    Matrix &X_, int x_s, int x_e, 
                    int y_s, int y_e) {
    //print_Matrix(X, x_s, x_e, y_s, y_e);
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            X[i][j] += X_[row][col];
        }
    }
}

void COPY_SUB_MATRIX(Matrix &X_,
                                   Matrix &X, int x_s,
                                   int x_e, int y_s, int y_e) {
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            X_[row][col] = X[i][j];
        }
    }
}

static inline void SUM_MATRIX(Matrix &Z,
                              Matrix &T) {
    for(int i = 0; i< Z.size(); ++i) {
        for(int j = 0; j < Z[0].size(); ++j) {
            Z[i][j] = Z[i][j] + T[i][j];
        }
    }
}

void Matrix_Multiply(Matrix &X, Matrix &Y, Matrix &Z, int x_row, int x_col, 
                     int y_row, int y_col, int z_row, int z_col, int n) {

    for(int i = 0; i<n; ++i){
        for(int j = 0; j<n; ++j) {
            for(int k = 0; k<n; ++k) {
                Z[z_row + i][z_col + j] += 
				    X[x_row + i][x_col + k] * Y[y_row + k][y_col + j];
            }
        }
    }
}


void init_work(Work *work, int size) {
    for(int i = 0; i< size; ++i){
        work[i].t_ = &PAR_REC_MEM;
        work[i].td_.X_ = &X;
        work[i].td_.Y_ = &Y;
        work[i].td_.Z_ = &Z;
	}
}

void set_work(Work *work, int x_row, int x_col, int y_row, int y_col, int z_row, 
	          int z_col, int n) {
	work->td_.x_row = x_row;
	work->td_.x_col = x_col;
	work->td_.y_row = y_row;
	work->td_.y_col = y_col;
	work->td_.z_row = z_row;
	work->td_.z_col = z_col;
    work->td_.n_ = n;
}

void PAR_REC_MEM(Matrix &X, Matrix &Y, Matrix &Z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, int id) {
        
    if(n == 1) {
        Matrix_Multiply(X, Y, Z, x_row, x_col, y_row, y_col, z_row, z_col, n);
    }
    else {
        Work *work = new Work[8];
		init_work(work, 8);
		set_work(&work[0], x_row, x_col, y_row, y_col, z_row, z_col, n/2);
		set_work(&work[1], x_row, x_col, y_row, y_col + n/2, z_row, z_col + n/2, n/2);
		set_work(&work[2], x_row + n/2, x_col, y_row, y_col, z_row + n/2, z_col, n/2);
		set_work(&work[3], x_row + n/2, x_col, y_row, y_col + n/2, z_row + n/2, z_col + n/2, n/2);
        pool[id]->push_back(work[0]);
        pool[id]->push_back(work[1]);
        pool[id]->push_back(work[2]);
        pool[id]->push_back(work[3]);

        set_work(&work[4], x_row, x_col + n/2, y_row + n/2, y_col, z_row, z_col, n/2);
		set_work(&work[5], x_row, x_col + n/2, y_row + n/2, y_col + n/2, z_row, z_col + n/2, n/2);
		set_work(&work[6], x_row + n/2, x_col + n/2, y_row + n/2, y_col, z_row + n/2, z_col, n/2);
		set_work(&work[7], x_row + n/2, x_col + n/2, y_row + n/2, y_col + n/2, z_row + n/2, z_col + n/2, n/2);
        pool[id]->push_back(work[4]);
        pool[id]->push_back(work[5]);
        pool[id]->push_back(work[6]);
        pool[id]->push_back(work[7]);

        /*
		PAR_REC_MEM(X, Y, Z, x_row, x_col, y_row, y_col, z_row, z_col, n/2);
        PAR_REC_MEM(X, Y, Z, x_row, x_col, y_row, y_col + n/2, z_row, z_col + n/2, n/2);
        PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col, y_row, y_col, z_row + n/2, z_col, n/2);
        PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col, y_row, y_col + n/2, z_row + n/2, z_col + n/2, n/2);

        PAR_REC_MEM(X, Y, Z, x_row, x_col + n/2, y_row + n/2, y_col, z_row, z_col, n/2);
        PAR_REC_MEM(X, Y, Z, x_row, x_col + n/2, y_row + n/2, y_col + n/2, z_row, z_col + n/2, n/2);
        PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col + n/2, y_row + n/2, y_col, z_row + n/2, z_col, n/2);
        PAR_REC_MEM(X, Y, Z, x_row + n/2, x_col + n/2, y_row + n/2, y_col + n/2, z_row + n/2, z_col + n/2, n/2);
		*/
    }
}

int main() {

    fillMatrix(X);
    fillMatrix(Y);

    // Print matrix
    printMatrix(X, 8);
    printMatrix(Y, 8);
    
    Work work;
    init_work(&work, 1);
    set_work(&work, 0, 0, 0, 0, 0, 0, 8);

    unsigned int num_threads = std::thread::hardware_concurrency();
	int i = 0;
    for(; i<num_threads; ++i) {
        pool.push_back(new Task(i));
    }
    
    // cilk spawn each thread in the pool
    for(int i = 0; i<num_threads; ++i) {
        pool[i]->push_back(work);
        pool[i]->run();
    }

    printMatrix(Z, 8);
    
	return 1;
}
