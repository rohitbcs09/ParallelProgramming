#include "mat_mul_sharing.h"
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

unsigned int num_threads = 1;

uint64_t g_seed = time(0);

// Thread Pool 
std::vector<Task *> pool;


void PAR_REC_MEM_BOTTOM_HALF(Matrix* X, Matrix* Y, Matrix* Z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, Sync* sync,
				 int id);
 
void PAR_REC_MEM(Matrix* X, Matrix* Y, Matrix* Z, int x_row, int x_col, 
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

void SUM_SUB_MATRIX(Matrix* X, 
                    Matrix* X_, int x_s, int x_e, 
                    int y_s, int y_e) {
    //print_Matrix(X, x_s, x_e, y_s, y_e);
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            (*X)[i][j] += (*X_)[row][col];
        }
    }
}

void COPY_SUB_MATRIX(Matrix* X_, Matrix* X, int x_s, int x_e, int y_s, int y_e) {
    for(int i = x_s, row = 0; i < x_e; ++i, ++row) {
        for(int j = y_s, col = 0; j < y_e; ++j, ++col) {
            (*X_)[row][col] = (*X)[i][j];
        }
    }
}

static inline void SUM_MATRIX(Matrix* Z, Matrix* T) {
    for(int i = 0; i< (*Z).size(); ++i) {
        for(int j = 0; j < Z[0].size(); ++j) {
            (*Z)[i][j] = (*Z)[i][j] + (*T)[i][j];
        }
    }
}


void Matrix_Multiply(Matrix *X, Matrix *Y, Matrix *Z, int x_row, int x_col, 
                    int y_row, int y_col, int z_row, int z_col, int n) {
    for(int i = 0; i<n; ++i){
        for(int j = 0; j<n; ++j) {
            for(int k = 0; k<n; ++k) {
                ((*Z)[z_row + i][z_col + j]) += 
                    ((*X)[x_row + i][x_col + k]) * ((*Y)[y_row + k][y_col + j]);
            }
        }
    }
}

void init_work(work *work, int end) {
    for(int i = 0; i<end; ++i){
        work[i].t_ = &PAR_REC_MEM;
        work[i].td_.X_ = &X;
        work[i].td_.Y_ = &Y;
        work[i].td_.Z_ = &Z;
        work[i].sync_  = NULL;
    }
}

void set_work(Work *work, int x_row, int x_col, int y_row, int y_col, int z_row, 
              int z_col, int n, Sync* sync) {
    work->td_.x_row = x_row;
    work->td_.x_col = x_col;
    work->td_.y_row = y_row;
    work->td_.y_col = y_col;
    work->td_.z_row = z_row;
    work->td_.z_col = z_col;
    work->td_.n_ = n;
	work->sync_ = sync;
}

Work* create_work(int x_row, int x_col, int y_row, int y_col, int z_row, 
                  int z_col , int n, Sync *sync) {
    Work* w = new Work;
    w->t_ = &PAR_REC_MEM;
    w->td_.X_ = &X;
    w->td_.Y_ = &Y;
    w->td_.Z_ = &Z;
	w->td_.x_row = x_row;
	w->td_.x_col = x_col;
	w->td_.y_row = y_row;
	w->td_.y_col = y_col;
	w->td_.z_row = z_row;
	w->td_.z_col = z_col;
	w->td_.n_ = n;
	w->sync_  = sync;
	return w;
}

void PAR_REC_MEM(Matrix *X, Matrix *Y, Matrix *Z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, Sync* sync,
                 int id) {
    std::cout << "ParRecMM " << " z_row: " << z_row << " z_col: " << z_col 
	          << " x_row: " << x_row << " x_col: " << x_col << " y_row: " << y_row 
			  << " y_col: " << y_col << " size: " << n;
	if(sync) {
	    std::cout << " syncType: " << sync->type_ << " syncValue: " << sync->val_<< "\n";
	}
	else {
	    std::cout << "\n";
    }
    if(n == 1) {
        Matrix_Multiply(X, Y, Z, x_row, x_col, y_row, y_col, z_row, z_col, 1);
		pool[id]->dec_sync_ref(sync);
		if(pool[id]->get_sync_ref_count(sync) == 0) {
		    //pool[id]->pop_sync();
			Sync* par = sync;
			while(par != NULL) {
		      	if(par->type_ == 1) {
				    if(pool[id]->get_sync_ref_count(par) == 0 ){
							Sync *new_sync = create_sync_state(par, 2, par->x_row, par->x_col,
								                               par->y_row, par->y_col, par->z_row,
															   par->z_col, par->n);
							pool[id]->push_sync(new_sync);
							PAR_REC_MEM_BOTTOM_HALF(X, Y, Z, par->x_row, par->x_col, 
		                                    par->y_row, par->y_col, par->z_row,
								        	par->z_col, par->n, new_sync, id);
							//pool[id]->dec_sync_ref(par);
							break;
				   	}
					else {
                        break;
					}
				}
		      	else {
				    if(pool[id]->get_sync_ref_count(par) == 0) {
						std::cout << "BEFORE PARENT :  " << par->val_ << " " << par <<"\n";
                        par = par->parent;
                        pool[id]->dec_sync_ref(par);
                        //if(pool[id]->get_sync_ref_count(par) <= 0) {
						//         pool[id]->pop_sync();
								 //continue;
					    std::cout << "AFTER PARENT :  " << par->val_ << " " << par << "\n";
                                 //break;
					    //}
					}
					else {
                        //pool[id]->dec_sync_ref(par);
                        break;
					}
				}
		    }
	    }
		return;
    }
    else {
		Sync* child_sync = create_sync_state(sync, 1, x_row, x_col, y_row, y_col, z_row, z_col, n);
        pool[id]->push_sync(child_sync);

        // pushing the work in Deque
		int rand_id =  fastrand() % num_threads;
        pool[rand_id]->push_back(create_work(x_row, x_col, y_row, y_col, z_row, 
		                                z_col, n/2, child_sync));

		rand_id =  fastrand() % num_threads;
        pool[rand_id]->push_back(create_work(x_row, x_col, y_row, y_col + n/2, 
                                        z_row, z_col + n/2, n/2, child_sync));

		rand_id =  fastrand() % num_threads;
        pool[id]->push_back(create_work(x_row + n/2, x_col, y_row, y_col, 
                                        z_row + n/2, z_col, n/2, child_sync));

        pool[id]->push_back(create_work(x_row + n/2, x_col, y_row, y_col + n/2, 
                                        z_row + n/2, z_col + n/2, n/2, child_sync));
    }
}

void PAR_REC_MEM_BOTTOM_HALF(Matrix *X, Matrix *Y, Matrix *Z, int x_row, int x_col, 
                 int y_row, int y_col, int z_row, int z_col, int n, Sync* sync,
                 int id) {
    std::cout << "Func2    " << " z_row: " << z_row << " z_col: " << z_col 
	          << " x_row: " << x_row << " x_col: " << x_col << " y_row: " << y_row 
			  << " y_col: " << y_col << " size: " << n;
	
	if(sync) {
	    std::cout << " syncType: " << sync->type_ << " syncValue: " << sync->val_<< "\n";
	}
	else {
	    std::cout << "\n";
    }

    // pushing the work in Deque
    int rand_id =  fastrand() % num_threads;
    pool[rand_id]->push_back(create_work(x_row, x_col + n/2, y_row + n/2, y_col, 
                             z_row, z_col, n/2, sync));
    rand_id =  fastrand() % num_threads;
    pool[rand_id]->push_back(create_work(x_row, x_col + n/2, y_row + n/2, y_col + n/2, 
                             z_row, z_col + n/2, n/2, sync));
    rand_id =  fastrand() % num_threads;
    pool[rand_id]->push_back(create_work(x_row + n/2, x_col + n/2, y_row + n/2, y_col,
                             z_row + n/2, z_col, n/2, sync));
    rand_id =  fastrand() % num_threads;
    pool[rand_id]->push_back(create_work(x_row + n/2, x_col + n/2, y_row + n/2, 
	                         y_col + n/2, z_row + n/2, z_col + n/2, n/2, sync));
}


int main() {

    fillMatrix(X);
    fillMatrix(Y);

    // Print matrix
    printMatrix(X, 8);
    printMatrix(Y, 8);

    num_threads = std::thread::hardware_concurrency();
    int i = 0;
    for(; i<num_threads; ++i) {
        pool.push_back(new Task(i));
    }
    
    pool[0]->push_back(create_work(0, 0, 0, 0, 0, 0, 8, NULL));

    // cilk spawn each thread in the pool
    for(int i = 0; i<num_threads; ++i) {
        pool[i]->run();
    }

    printMatrix(Z, 8);
    
    return 1;
}

