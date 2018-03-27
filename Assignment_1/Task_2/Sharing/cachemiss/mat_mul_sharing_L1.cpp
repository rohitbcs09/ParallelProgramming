#include "mat_mul_sharing.h"
#include <ctime>
#include <stdint.h> 
#include <cilk/cilk.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <papi.h>


using namespace std;

/*

FAST RANDOM NUMBER GENERATOR:
============================

https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c

*/

void handle_error(int err){
    std::cerr << "PAPI error: " << err << std::endl;
}

unsigned int cores = 1;

uint64_t g_seed = time(0);

// Thread Pool 
std::vector<Task *> pool;


void PAR_REC_MEM_BOTTOM_HALF(Matrix* x, Matrix* y, Matrix* z, int x_row, 
                             int x_col, int y_row, int y_col, int z_row, 
                             int z_col, int n, Sync* sync, int id);
 
void PAR_REC_MEM(Matrix* x, Matrix* y, Matrix* z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, Sync* sync,
                 int id);
 
void fillMatrix(Matrix &arr, int n) {
    for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            arr[i][j] = fastrand();
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
        if(cur && cur->s_type == 1 && pool[id]->get_sync_ref_count(cur) == 0) {
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
        else if(cur && cur->s_type == 2 && pool[id]->get_sync_ref_count(cur) == 0) {
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
    if(n == 16) {
        // Matrix_Multiply(x, y, z, x_row, x_col, y_row, y_col, z_row, z_col, 1);
        for(int i = 0; i<n; ++i){
            for(int k = 0; k<n; ++k) {
                for(int j = 0; j<n; ++j) {
                    ((*z)[z_row + i][z_col + j]) += 
                        ((*x)[x_row + i][x_col + k]) * ((*y)[y_row + k][y_col + j]);
                }
            }
        }
        pool[id]->dec_sync_ref_count(sync);
        if(sync && pool[id]->get_sync_ref_count(sync) == 0) {
            pool[id]->pop_sync();
            //update_sync_queue_bottom_half(x, y, z, sync, id);
            Sync* cur = sync;
            while(cur != NULL) {
                if(cur && cur->s_type == 1 && pool[id]->get_sync_ref_count(cur) == 0) {
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
                else if(cur && cur->s_type == 2 && pool[id]->get_sync_ref_count(cur) == 0) {
                    cur = cur->s_parent;
                    pool[id]->dec_sync_ref_count(cur);
                }
                else {
                    break;
                }
            }
        }
        return;
    }
    else {
        Sync* child_sync = create_sync_state(sync, 1, x_row, x_col, y_row, 
                                             y_col, z_row, z_col, n);
        pool[id]->push_sync(child_sync);

        // pushing the work in Deque
        int rand_id =  fastrand() % cores;
        pool[rand_id]->push_front(create_work(x, y, z, x_row, x_col, y_row, 
                                             y_col, z_row, z_col, n/2, 
                                             child_sync));
        rand_id =  fastrand() % cores;
        pool[rand_id]->push_front(create_work(x, y, z, x_row, x_col, y_row, y_col + 
                                        n/2, z_row, z_col + n/2, n/2, 
                                        child_sync));
        rand_id =  fastrand() % cores;
        pool[rand_id]->push_front(create_work(x, y, z, x_row + n/2, x_col, y_row, 
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
    int rand_id =  fastrand() % cores;
    pool[rand_id]->push_front(create_work(x, y, z, x_row, x_col + n/2, y_row + 
                                         n/2, y_col, z_row, z_col, n/2, sync));
    rand_id =  fastrand() % cores;
    pool[rand_id]->push_front(create_work(x, y, z, x_row, x_col + n/2, y_row + 
                                         n/2, y_col + n/2, z_row, z_col + n/2, 
                                         n/2, sync));
    rand_id =  fastrand() % cores;
    pool[rand_id]->push_front(create_work(x, y, z, x_row + n/2, x_col + n/2, 
                                         y_row + n/2, y_col, z_row + n/2, 
                                         z_col, n/2, sync));
    rand_id =  fastrand() % cores;
    pool[rand_id]->push_back(create_work(x, y, z, x_row + n/2, x_col + n/2, 
                                         y_row + n/2, y_col + n/2, z_row + n/2,
                                         z_col + n/2, n/2, sync));
}

int main(int argc, char* argv[]) {
    if(argc < 2)  {
        std::cout << "Please peovide the input matrix size as 1st arg\n";
        return 1;
    }
    
    int n = atoi(argv[1]);
    
    if(argc == 3) {
       cores =  atoi(argv[2]);
    }

    // Input Matrices
    Matrix X(n, std::vector<uint64_t>(n, 0));
    Matrix Y(n, std::vector<uint64_t>(n, 0));
    Matrix Z(n, std::vector<uint64_t>(n, 0));

    fillMatrix(X, n);
    fillMatrix(Y, n);

    // Print matrix
    //printMatrix(X, n);
    //printMatrix(Y, n);

    //num_threads = std::thread::hardware_concurrency();
    if(n < 32) {
        //printMatrix(X, n);
        //printMatrix(Y, n);
        Matrix_Multiply(&X, &Y, &Z, 0, 0, 0, 0, 0, 0, n);
        //printMatrix(Z, n);
        return 1;
    }
    int i = 0;
    for(; i<cores; ++i) {
        pool.push_back(new Task(i, cores));
    }
    pool[0]->push_back(create_work(&X, &Y, &Z, 0, 0, 0, 0, 0, 0, n, NULL));
    using namespace std::chrono;    
    high_resolution_clock::time_point start_time = high_resolution_clock::now();

    int numEvents = 1;
    long long values[1];
    int events[1] = {PAPI_L1_TCM};

    if (PAPI_start_counters(events, numEvents) != PAPI_OK) {
            handle_error(1);
    }

    // cilk spawn each thread in the pool
    for(int i = 0; i<cores-1; ++i) {
        cilk_spawn pool[i]->run();
    }
    pool[0]->run();
    cilk_sync;
    high_resolution_clock::time_point end_time = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end_time - start_time);

    if ( PAPI_stop_counters(values, numEvents) != PAPI_OK) {
            handle_error(1);
    }

    std::cout<<"L1 misses: "<<values[0] << " for size " <<  argv[1] << " and cores " << argv[2] <<  std::endl;

    // std::cout << "Exectution Time: " << time_span.count() << " seconds.";
    std::cout << std::endl;

    //printMatrix(Z, n);
    
    return 1;
}

