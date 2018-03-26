#include "mat_mul_steal.h"
#include <ctime>
#include <stdint.h>
using namespace std;

/*
FAST RANDOM NUMBER GENERATOR:
============================

https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c

*/

#define ROWS 8 
#define COLS 8 

unsigned int num_threads = 1;

uint64_t g_seed = time(0);

// Thread Pool 
std::vector<Task *> pool;


void PAR_REC_MEM_BOTTOM_HALF(Matrix* x, Matrix* y, Matrix* z, int x_row, 
                             int x_col, int y_row, int y_col, int z_row, 
							 int z_col, int n, Sync* sync, int id);
 
void PAR_REC_MEM(Matrix* x, Matrix* y, Matrix* z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, Sync* sync,
                 int id);
 
void fillMatrix(Matrix &arr) {
    for(int i = 0; i<ROWS; ++i) {
        for(int j = 0; j<COLS; ++j) {
            arr[i][j] = j+1; //fastrand() % 10;
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

void SUM_SUB_MATRIX(Matrix* x, 
                    Matrix* x_, int x_s, int x_e, 
                    int y_s, int y_e) {
    //print_Matrix(x, x_s, x_e, y_s, y_e);
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            (*x)[i][j] += (*x_)[row][col];
        }
    }
}

void COPY_SUB_MATRIX(Matrix* x_, Matrix* x, int x_s, int x_e, int y_s, 
                     int y_e) {
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            (*x_)[row][col] = (*x)[i][j];
        }
    }
}

static inline void SUM_MATRIX(Matrix* z, Matrix* T) {
    for(int i = 0; i< (*z).size(); ++i) {
        for(int j = 0; j < z[0].size(); ++j) {
            (*z)[i][j] = (*z)[i][j] + (*T)[i][j];
        }
    }
}


void Matrix_Multiply(Matrix *x, Matrix *y, Matrix *z, int x_row, int x_col, 
                    int y_row, int y_col, int z_row, int z_col, int n) {
    for(int i = 0; i<n; ++i){
      for(int k = 0; k<n; ++k) {
        for(int j = 0; j<n; ++j) {
                ((*z)[z_row + i][z_col + j]) += 
                    ((*x)[x_row + i][x_col + k]) * ((*y)[y_row + k][y_col + j]);
            }
        }
    }
}

Work* create_work(Matrix* x, Matrix* y, Matrix* z, int x_row, int x_col,
                  int y_row, int y_col, int z_row, int z_col , int n, 
                  Sync *sync) {
    Work* w = new Work;
    w->w_t = &PAR_REC_MEM;
    w->w_td.t_X = x;
    w->w_td.t_Y = y;
    w->w_td.t_Z = z;
    w->w_td.t_x_row = x_row;
    w->w_td.t_x_col = x_col;
    w->w_td.t_y_row = y_row;
    w->w_td.t_y_col = y_col;
    w->w_td.t_z_row = z_row;
    w->w_td.t_z_col = z_col;
    w->w_td.t_n = n;
    w->w_sync  = sync;
    return w;
}

void update_sync_queue_bottom_half(Matrix *x, Matrix *y, Matrix *z, Sync *cur, 
                                   int id) {
    while(cur != NULL) {
        if(cur->s_type == 1 && pool[id]->get_sync_ref_count(cur) == 0) {
            Sync *sync_2 = create_sync_state(cur->s_parent, 2, cur->s_x_row, 
                                             cur->s_x_col, cur->s_y_row, 
                                             cur->s_y_col, cur->s_z_row,
                                             cur->s_z_col, cur->s_n);
            pool[id]->push_sync(sync_2);
            PAR_REC_MEM_BOTTOM_HALF(x, y, z, sync_2->s_x_row, sync_2->s_x_col, 
                                    sync_2->s_y_row, sync_2->s_y_col, 
                                    sync_2->s_z_row, sync_2->s_z_col, 
                                    sync_2->s_n, sync_2, id);
            break;
        }
        else if(cur->s_type == 2 && pool[id]->get_sync_ref_count(cur) == 0) {
            cur = cur->s_parent;
            pool[id]->dec_sync_ref_count(cur);
        }
        else {
            break;
        }
    }
}

void PAR_REC_MEM(Matrix *x, Matrix *y, Matrix *z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, Sync* sync,
                 int id) {
    if(n == 1) {
        Matrix_Multiply(x, y, z, x_row, x_col, y_row, y_col, z_row, z_col, 1);
        pool[id]->dec_sync_ref_count(sync);
        if(pool[id]->get_sync_ref_count(sync) == 0) {
            pool[id]->pop_sync();
            update_sync_queue_bottom_half(x, y, z, sync, id);
        }
        return;
    }
    else {
        Sync* child_sync = create_sync_state(sync, 1, x_row, x_col, y_row, 
                                             y_col, z_row, z_col, n);
        pool[id]->push_sync(child_sync);

        // pushing the work in Deque
        pool[id]->push_back(create_work(x, y, z, x_row, x_col, y_row, y_col, 
                                        z_row, z_col, n/2, child_sync));
        pool[id]->push_back(create_work(x, y, z, x_row, x_col, y_row, y_col + 
                                        n/2, z_row, z_col + n/2, n/2, 
                                        child_sync));
        pool[id]->push_back(create_work(x, y, z, x_row + n/2, x_col, y_row, 
                                        y_col, z_row + n/2, z_col, n/2, 
                                        child_sync));
        pool[id]->push_back(create_work(x, y, z, x_row + n/2, x_col, y_row, 
                                        y_col + n/2, z_row + n/2, z_col + n/2, 
                                        n/2, child_sync));
    }
    return;
}

void PAR_REC_MEM_BOTTOM_HALF(Matrix *x, Matrix *y, Matrix *z, int x_row, 
                             int x_col, int y_row, int y_col, int z_row, 
                             int z_col, int n, Sync* sync, int id) {

    // pushing the work in Deque
    pool[id]->push_back(create_work(x, y, z, x_row, x_col + n/2, y_row + n/2, 
                                    y_col, z_row, z_col, n/2, sync));
    pool[id]->push_back(create_work(x, y, z, x_row, x_col + n/2, y_row + n/2, 
                                    y_col + n/2, z_row, z_col + n/2, n/2, 
                                    sync));
    pool[id]->push_back(create_work(x, y, z, x_row + n/2, x_col + n/2, y_row + 
                                    n/2, y_col, z_row + n/2, z_col, n/2, 
                                    sync));
    pool[id]->push_back(create_work(x, y, z, x_row + n/2, x_col + n/2, y_row + 
                                    n/2, y_col + n/2, z_row + n/2, z_col + n/2, 
                                    n/2, sync));
}

Work* Task::steal_random_work() {
    int id = fastrand() % num_threads;
    return pool[id]->pop_front();
}


int main(int argc, char* argv[]) {
    if(argc < 2)  {
        std::cout << "Please peovide the input matrix size as 1st arg\n";
		return 1;
	}
	
	int n = atoi(argv[1]);

    // Input Matrices
    Matrix X(n, std::vector<uint64_t>(n, 0));
    Matrix Y(n, std::vector<uint64_t>(n, 0));
    Matrix Z(n, std::vector<uint64_t>(n, 0));

    fillMatrix(X);
    fillMatrix(Y);

    // Print matrix
    printMatrix(X, n);
    printMatrix(Y, n);

    num_threads = std::thread::hardware_concurrency();

    int i = 0;
    for(; i<num_threads; ++i) {
        pool.push_back(new Task(i, num_threads));
    }
    
    pool[0]->push_back(create_work(&X, &Y, &Z, 0, 0, 0, 0, 0, 0, n, NULL));

    // cilk spawn each thread in the pool
    for(int i = 0; i<num_threads; ++i) {
        pool[i]->run();
    }

    printMatrix(Z, n);
    
    return 1;
}

