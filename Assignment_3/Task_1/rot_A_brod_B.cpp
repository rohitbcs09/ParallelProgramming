#include <iostream>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <chrono>

using namespace std;

// Input Matrices
int **X; 
int **Y;
int **Z;

int g_seed;

int processors = 1;

int fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
}

void create_2d_array(int ***array, int n, int m) {

    int *p = (int *)malloc(n*m*sizeof(int));
    if(!p) {
        std::cout << "Malloc Failed !!!\n";
        return;
    }

    (*array) = (int**)malloc(n*sizeof(int*));
    if (!(*array)) {
       free(p);
       std::cout << "Malloc Failed !!!\n";
       return;
    }

    for (int i=0; i<n; i++) {
       (*array)[i] = &(p[i*m]);
    }
    return;
}

void delete_2d_array(int ***array) {
    free(&((*array)[0][0]));
    free(*array);
    return;
}

void fillMatrix(int **arr, int n) {
    int count = 0;
    for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            arr[i][j] = ++count; //fastrand();
        }
    }
}

void fillMatrixVal(int **arr, int n, int val) {
    for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            arr[i][j] = val; //fastrand();
        }
    }
}

void printMatrix(int **arr, int n) {
     std::cout << "\n";
     for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }   
}

void init_sub_matrix(int **arr, int m, int n) {
    for(int i = 0; i< m; ++i) {
        for(int j = 0; j< m; ++j) {
            arr[i][j] = 0;
        }
    }
}


void Matrix_Multiply(int **X, int **Y, int **Z, int x_row, int x_col, 
                     int y_row, int y_col, int z_row, int z_col, int n) {
    for(int i = 0; i<n; ++i){
        for(int k = 0; k<n; ++k) {
            for(int j = 0; j<n; ++j) {
                Z[z_row + i][z_col + j] += X[x_row + i][x_col + k] * Y[y_row + k][y_col + j];
            }
        }
    }
}

void MM_rotate_A_broadcast_B(int n, int p, int rank) {
    int *mat_X = NULL;
    int *mat_Y = NULL;
    int *mat_Z = NULL;
    int *g_copy_mat = NULL;

    int sqrt_p = std::sqrt(p);
    int proc_mat_size = n / sqrt_p;
    int mat_size[2] = {n , n};
    int sub_mat_size[2] = {proc_mat_size, proc_mat_size};
    int mat_start[2] = {0, 0};

    int **sub_matrix_X;
    create_2d_array(&sub_matrix_X, proc_mat_size, proc_mat_size);

    int **sub_matrix_Y;
    create_2d_array(&sub_matrix_Y, proc_mat_size, proc_mat_size);

    int **sub_matrix_Z;
    create_2d_array(&sub_matrix_Z, proc_mat_size, proc_mat_size);
    init_sub_matrix(sub_matrix_Z, proc_mat_size, proc_mat_size);

    MPI_Status status;
    MPI_Datatype mat_type, sub_mat_type;
    MPI_Type_create_subarray(2, mat_size, sub_mat_size, mat_start, MPI_ORDER_C, MPI_INT, &mat_type);
    MPI_Type_create_resized(mat_type, 0, proc_mat_size * sizeof(int), &sub_mat_type);
    MPI_Type_commit(&sub_mat_type);
    
    int send_mat_count[p]; 
    int send_mat_start[p];

    if(rank == 0) {
        int start = 0;
        for(int i=0; i<sqrt_p; ++i){
            for(int j=0; j<sqrt_p; ++j){
                send_mat_count[start]=1;
                send_mat_start[start++]= (sqrt_p * proc_mat_size * i) + j;
                //std::cout << (sqrt_p * proc_mat_size * i) + j << " ";
            }
        }
        /*for(int i = 0; i< sqrt_p * sqrt_p; ++i) {
            std::cout << send_mat_start[i] << " ";
        }
        std::cout << "\n";*/
        mat_X = &(X[0][0]);
        mat_Y = &(Y[0][0]);
        mat_Z = &(Z[0][0]);
    }

    MPI_Scatterv(mat_X, send_mat_count, send_mat_start, sub_mat_type, &(sub_matrix_X[0][0]),
                 p, MPI_INT, 0, MPI_COMM_WORLD);
 
    MPI_Scatterv(mat_Y, send_mat_count, send_mat_start, sub_mat_type, &(sub_matrix_Y[0][0]),
                 p, MPI_INT, 0, MPI_COMM_WORLD);
 
    MPI_Scatterv(mat_Z, send_mat_count, send_mat_start, sub_mat_type, &(sub_matrix_Z[0][0]),
                 p, MPI_INT, 0, MPI_COMM_WORLD);
          
    /*std::cout << "Rank : " << rank << "\n";
    for(int j = 0; j< proc_mat_size; ++j) {
        for(int k =0; k< proc_mat_size; ++k) {
            std::cout << sub_matrix_X[j][k] << " ";
        }
        std::cout << "\n";
    }*/

    int color = rank % sqrt_p;

    MPI_Comm col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &col_comm);


    MPI_Barrier(MPI_COMM_WORLD);

    int tag = 99;
    int p_src_row = rank / sqrt_p; 
    int p_src_col = rank % sqrt_p;

    int x_send_to_p_dst_row = p_src_row;     
    int x_recv_from_p_src_row = p_src_row;

    int **bcast_ptr;
    create_2d_array(&bcast_ptr, proc_mat_size, proc_mat_size);
    init_sub_matrix(bcast_ptr, proc_mat_size, proc_mat_size);
    int *g_bcast_ptr = &(bcast_ptr[0][0]);
	
    for(int l = 0; l < sqrt_p; l++) {
        
        // Broad cast
        int k = l + p_src_col;
        //std::cout << "j = " << p_src_col << " l = " << l << " i = " << p_src_row << " bcast rank = "  << (rank - p_src_col) / n << std::endl;
        if (k == p_src_row) {
            g_bcast_ptr = &(sub_matrix_Y[0][0]);;
        }
	
	MPI_Bcast(g_bcast_ptr, proc_mat_size * proc_mat_size, MPI_INT, p_src_row, col_comm);

        Matrix_Multiply(sub_matrix_X, &g_bcast_ptr, sub_matrix_Z,
                        0, 0, 0, 0, 0, 0, proc_mat_size);

        int x_send_to_p_dst_col = (sqrt_p + p_src_col - 1) % sqrt_p; 
        int x_p_dst_rank = (sqrt_p * x_send_to_p_dst_row ) + x_send_to_p_dst_col;

        int x_recv_from_p_src_col = (sqrt_p + p_src_col + 1) % sqrt_p; 
        int x_p_src_rank = (sqrt_p * x_recv_from_p_src_row ) + x_recv_from_p_src_col;

        MPI_Sendrecv_replace(&(sub_matrix_X[0][0]), proc_mat_size * proc_mat_size, MPI_INT,
                                  x_p_dst_rank, tag, x_p_src_rank, tag, MPI_COMM_WORLD, &status); 

    }

    std::cout << "Hello" << std::endl;

    MPI_Gatherv(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size,  MPI_INT, mat_Z, send_mat_count, 
                send_mat_start, sub_mat_type, 0, MPI_COMM_WORLD);
    if(rank == 0) {
        std::cout << "\nMASTER OUTPUT:\n";
        printMatrix(Z, n);
    }

    return;
}

int main(int argc, char *argv[]) {
    if(argc < 2)  {
        std::cout << "Missing Input params - ibrun -n <processors> <binary> <matrix_size>\n";
        return 1;
    }
    int n = atoi(argv[1]);

    srand(time(NULL));
    g_seed=rand();

    int myrank, v = 121;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);
    
    if(myrank == 0) {
        create_2d_array(&X, n, n);
        fillMatrix(X, n);
        printMatrix(X, n);

        create_2d_array(&Y, n, n);
        fillMatrix(Y, n);

        create_2d_array(&Z, n, n);
        init_sub_matrix(Z, n, n);
    }


    using namespace std::chrono;    
    high_resolution_clock::time_point start_time = high_resolution_clock::now();

    MM_rotate_A_broadcast_B(n, processors, myrank);

    high_resolution_clock::time_point end_time = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end_time - start_time);

    if(myrank == 0) {
        std::cout << "Exectution Time: " << time_span.count() << " seconds.";
        std::cout << std::endl;
    }
    return 1;
}
