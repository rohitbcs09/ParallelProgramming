#include <iostream>
#include <vector>
#include <ctime>
#include <stdint.h>
using namespace std;

/*

FAST RANDOM NUMBER GENERATOR:
    o https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c

*/

#define ROWS 4 
#define COLS 4 

uint64_t g_seed = time(0);

static inline uint64_t fastrand() { 
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
} 

void fillMatrix(std::vector<std::vector<uint64_t> > &arr) {
    for(int i = 0; i<ROWS; ++i) {
        for(int j = 0; j<COLS; ++j) {
            arr[i][j] = fastrand() % 10;
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
                     std::vector<std::vector<uint64_t> > &Z, int n) {

    for(int i = 0; i<n; ++i){
        for(int j = 0; j<n; ++j) {
            Z[i][j] = 0;
            for(int k = 0; k<n; ++k) {
                Z[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
}

std::vector<std::vector<uint64_t> >
    PAR_REC_MEM(std::vector<std::vector<uint64_t> > X,
                std::vector<std::vector<uint64_t> > Y, int n) {
        
    std::vector<std::vector<uint64_t> > Z(n, std::vector<uint64_t>(n, 0));
    if(n == 1) {
        Matrix_Multiply(X, Y, Z, n);
    }
    else {
        std::vector<std::vector<uint64_t> > X_11(n/2, std::vector<uint64_t>(n/2,0));
        std::vector<std::vector<uint64_t> > X_12(n/2, std::vector<uint64_t>(n/2,0));
        std::vector<std::vector<uint64_t> > X_21(n/2, std::vector<uint64_t>(n/2,0));
        std::vector<std::vector<uint64_t> > X_22(n/2, std::vector<uint64_t>(n/2,0));


        std::vector<std::vector<uint64_t> > Y_11(n/2, std::vector<uint64_t>(n/2,0));
        std::vector<std::vector<uint64_t> > Y_12(n/2, std::vector<uint64_t>(n/2,0));
        std::vector<std::vector<uint64_t> > Y_21(n/2, std::vector<uint64_t>(n/2,0));
        std::vector<std::vector<uint64_t> > Y_22(n/2, std::vector<uint64_t>(n/2,0));

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

        std::vector<std::vector<uint64_t> > Z_11 = PAR_REC_MEM(X_11, Y_11, n/2);
        std::vector<std::vector<uint64_t> > T_11 = PAR_REC_MEM(X_12, Y_21, n/2);

        std::vector<std::vector<uint64_t> > Z_12 = PAR_REC_MEM(X_11, Y_12, n/2);
        std::vector<std::vector<uint64_t> > T_12 = PAR_REC_MEM(X_12, Y_22, n/2);

        std::vector<std::vector<uint64_t> > Z_21 = PAR_REC_MEM(X_21, Y_11, n/2);
        std::vector<std::vector<uint64_t> > T_21 = PAR_REC_MEM(X_22, Y_21, n/2);

        std::vector<std::vector<uint64_t> > Z_22 = PAR_REC_MEM(X_21, Y_12, n/2);
        std::vector<std::vector<uint64_t> > T_22 = PAR_REC_MEM(X_22, Y_22, n/2);

        SUM_SUB_MATRIX(Z, Z_11, 0, n/2, 0, n/2);
        SUM_SUB_MATRIX(Z, Z_12, 0, n/2, n/2, n);
        SUM_SUB_MATRIX(Z, Z_21, n/2, n, 0, n/2);
        SUM_SUB_MATRIX(Z, Z_22, n/2, n, n/2, n);

        SUM_SUB_MATRIX(Z, T_11, 0, n/2, 0, n/2);
        SUM_SUB_MATRIX(Z, T_12, 0, n/2, n/2, n);
        SUM_SUB_MATRIX(Z, T_21, n/2, n, 0, n/2);
        SUM_SUB_MATRIX(Z, T_22, n/2, n, n/2, n);
    }
    return Z;
}

int main() {
    std::vector< std::vector<uint64_t> > A(ROWS, std::vector<uint64_t>(COLS, 0));
    std::vector< std::vector<uint64_t> > B(ROWS, std::vector<uint64_t>(COLS, 0));
    
    fillMatrix(A);
    fillMatrix(B);

    // Print array
    printMatrix(A, 4);
    printMatrix(B, 4);

    std::vector< std::vector<uint64_t> > Z = PAR_REC_MEM(A, B, 4);
    printMatrix(Z, 4);

    return 1;
}

