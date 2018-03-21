#include "../include/rec_mat_mul.h"
#include <ctime>
#include <stdint.h>
using namespace std;

/*

FAST RANDOM NUMBER GENERATOR:
    o https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c

*/

#define ROWS 4 
#define COLS 4 

// Input Matrices
Matrix A(ROWS, std::vector<uint64_t>(COLS, 0));
Matrix B(ROWS, std::vector<uint64_t>(COLS, 0));
Matrix T(ROWS, std::vector<uint64_t>(COLS, 0));

// Thread Pool 
std::vector<Task *> pool;

uint64_t g_seed = time(0);

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

void Matrix_Multiply(Matrix &X, Matrix &Y, Matrix &Z, int n) {
    for(int i = 0; i<n; ++i){
        for(int j = 0; j<n; ++j) {
            Z[i][j] = 0;
            for(int k = 0; k<n; ++k) {
                Z[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
}

void PAR_REC_MEM(Matrix X, Matrix Y, Matrix Z, int n, int id) {
        
    if(n == 1) {
        //Matrix_Multiply(X, Y, Z, n);
	    std::cout << "Base Case " << "\n";

    }
    else {
	    std::cout << "Inside Rec function " << "\n";
        Matrix X_11(n/2, std::vector<uint64_t>(n/2,0));
        Matrix X_12(n/2, std::vector<uint64_t>(n/2,0));
        Matrix X_21(n/2, std::vector<uint64_t>(n/2,0));
        Matrix X_22(n/2, std::vector<uint64_t>(n/2,0));

        Matrix Y_11(n/2, std::vector<uint64_t>(n/2,0));
        Matrix Y_12(n/2, std::vector<uint64_t>(n/2,0));
        Matrix Y_21(n/2, std::vector<uint64_t>(n/2,0));
        Matrix Y_22(n/2, std::vector<uint64_t>(n/2,0));

        Matrix Z_11(n/2, std::vector<uint64_t>(n/2,0));
        Matrix Z_12(n/2, std::vector<uint64_t>(n/2,0));
        Matrix Z_21(n/2, std::vector<uint64_t>(n/2,0));
        Matrix Z_22(n/2, std::vector<uint64_t>(n/2,0));

        Matrix T_11(n/2, std::vector<uint64_t>(n/2,0));
        Matrix T_12(n/2, std::vector<uint64_t>(n/2,0));
        Matrix T_21(n/2, std::vector<uint64_t>(n/2,0));
        Matrix T_22(n/2, std::vector<uint64_t>(n/2,0));

        // Creating sub-matrix of X
        COPY_SUB_MATRIX(X_11, X, 0, n/2, 0, n/2);
        COPY_SUB_MATRIX(X_12, X, 0, n/2, n/2, n);
        COPY_SUB_MATRIX(X_21, X, n/2, n, 0, n/2);
        COPY_SUB_MATRIX(X_22, X, n/2, n, n/2, n);
        
        // Creating sub-matrix of Y 
        COPY_SUB_MATRIX(Y_11, Y, 0, n/2, 0, n/2);
        COPY_SUB_MATRIX(Y_12, Y, 0, n/2, n/2, n);
        COPY_SUB_MATRIX(Y_21, Y, n/2, n, 0, n/2);
        COPY_SUB_MATRIX(Y_22, Y, n/2, n, n/2, n);
       
		/*Work top_work[4];

		top_work[0].set_work(&PAR_REC_MEM, &X_11, &Y_11, &Z_11, n/2);
		top_work[1].set_work(&PAR_REC_MEM, &X_12, &Y_12, &Z_12, n/2);
		top_work[2].set_work(&PAR_REC_MEM, &X_21, &Y_21, &Z_21, n/2);
		top_work[3].set_work(&PAR_REC_MEM, &X_21, &Y_12, &Z_22, n/2);*/

        Work work;
        work.t_ = &PAR_REC_MEM;
        work.td_.X_ = &A;
        work.td_.Y_ = &B;
        work.td_.Z_ = &T;
        work.td_.n_ = n/2;
        pool[id]->push_back(work);
        pool[id]->push_back(work);
        pool[id]->push_back(work);
        pool[id]->push_back(work);

        //PAR_REC_MEM(X_21, Y_12, Z_22, n/2, id);

        // sync 1 

		/*Work bottom_work[4];
		bottom_work[0].set_work(&PAR_REC_MEM, &X_12, &Y_21, &T_11, n/2);
		bottom_work[1].set_work(&PAR_REC_MEM, &X_12, &Y_22, &T_12, n/2);
		bottom_work[2].set_work(&PAR_REC_MEM, &X_22, &Y_21, &T_21, n/2);
		bottom_work[3].set_work(&PAR_REC_MEM, &X_22, &Y_22, &T_22, n/2);*/
        pool[id]->push_back(work);
        pool[id]->push_back(work);
        pool[id]->push_back(work);
        pool[id]->push_back(work);

        //PAR_REC_MEM(X_22, Y_22, T_22, n/2, id);

        // sync 2

        SUM_SUB_MATRIX(Z, Z_11, 0, n/2, 0, n/2);
        SUM_SUB_MATRIX(Z, Z_12, 0, n/2, n/2, n);
        SUM_SUB_MATRIX(Z, Z_21, n/2, n, 0, n/2);
        SUM_SUB_MATRIX(Z, Z_22, n/2, n, n/2, n);

        SUM_SUB_MATRIX(Z, T_11, 0, n/2, 0, n/2);
        SUM_SUB_MATRIX(Z, T_12, 0, n/2, n/2, n);
        SUM_SUB_MATRIX(Z, T_21, n/2, n, 0, n/2);
        SUM_SUB_MATRIX(Z, T_22, n/2, n, n/2, n);
    }
}

int main() {

    fillMatrix(A);
    fillMatrix(B);

    // Print matrix
    printMatrix(A, 4);
    printMatrix(B, 4);
    
    Work work;
    work.t_ = &PAR_REC_MEM;
    work.td_.X_ = &A;
    work.td_.Y_ = &B;
    work.td_.Z_ = &T;
    work.td_.n_ = 4;

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

    printMatrix(T, 4);
    
	return 1;
}
