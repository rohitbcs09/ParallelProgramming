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
            arr[i][j] = /*++count; */ fastrand();
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

void update_result(int *temp, int rank, int **Z, int n , int p) {
    int sqrt_p = std::sqrt(p);
    int proc_mat_size = n / sqrt_p;
 
    int row = rank / sqrt_p;
    int col = rank % sqrt_p;
    int *trav = temp;

    for (int i = row * proc_mat_size; i < row * proc_mat_size + proc_mat_size; ++i) {
        for (int j = col * proc_mat_size; j < col * proc_mat_size + proc_mat_size; ++j) {

            //std::cout << " i = " << i << " j = " << j << " val = " << *trav;
            Z[i][j] += *trav;
            trav++;
        }
        //std::cout << "\n";
    }
}

int * get_buff_copy(int rank,int  **arr, int size, int procs) {
    int sqrt_p = std::sqrt(procs);
    int proc_mat_size = size / sqrt_p;

    int * res = (int *)malloc(proc_mat_size*proc_mat_size*sizeof(int));
    int *trav = res;
    int row = rank / sqrt_p;
    int col = rank % sqrt_p;

    for (int i = row * proc_mat_size; i < row * proc_mat_size + proc_mat_size; ++i) {
        for (int j = col * proc_mat_size; j < col * proc_mat_size + proc_mat_size; ++j) {
            *trav = arr[i][j];
            trav++;
        }
    }

    return res;
}

void MM_rotate_A_rotate_B(int n, int p, int rank) {
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
    

    MPI_Status x_status[n];
    MPI_Status y_status[n];
    MPI_Status z_status[n];

    int x_tag = 40;
    int y_tag = 41;
    int z_tag = 42;

    if(rank == 0) {

	for (int i = 0; i < proc_mat_size; ++i) {
	   for (int j = 0; j < proc_mat_size; ++j) {
               sub_matrix_X[i][j] = *(*X + i * n + j);
	       sub_matrix_Y[i][j] = *(*Y + i * n + j);
	       sub_matrix_Z[i][j] = *(*Z + i * proc_mat_size + j);
	   }
	}


        // send to all processors
        for (int i = 1; i < sqrt_p * sqrt_p; ++i) {
            int *sub_X = get_buff_copy(i, X, n, p);
            //MPI_Isend(sub_X, proc_mat_size * proc_mat_size, MPI_INT, i, x_tag, MPI_COMM_WORLD, &x_request[i]);
            MPI_Send(sub_X, proc_mat_size * proc_mat_size, MPI_INT, i, x_tag, MPI_COMM_WORLD);

            int *sub_Y = get_buff_copy(i, Y, n, p);
            //MPI_Isend(sub_Y, proc_mat_size * proc_mat_size, MPI_INT, i, y_tag, MPI_COMM_WORLD, &y_request[i]);
            MPI_Send(sub_Y, proc_mat_size * proc_mat_size, MPI_INT, i, y_tag, MPI_COMM_WORLD);

            int *sub_Z = get_buff_copy(i, Z, n, p);
            //MPI_Isend(sub_Z, proc_mat_size * proc_mat_size, MPI_INT, i, z_tag, MPI_COMM_WORLD, &z_request[i]);
            MPI_Send(sub_Z, proc_mat_size * proc_mat_size, MPI_INT, i, z_tag, MPI_COMM_WORLD);

            free(sub_X);
            free(sub_Y);
            free(sub_Z);
        }

    } 

    if (rank != 0) {
        //MPI_Irecv(&(sub_matrix_X[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_request[rank]);
        MPI_Recv(&(sub_matrix_X[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_status[rank]);
        //MPI_Irecv(&(sub_matrix_Y[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_request[rank]);
        MPI_Recv(&(sub_matrix_Y[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_status[rank]);
        //MPI_Irecv(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_request[rank]);
        MPI_Recv(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_status[rank]);
    
        //MPI_Wait(&x_request[rank], &x_status[rank]);
        //MPI_Wait(&y_request[rank], &y_status[rank]);
        //MPI_Wait(&z_request[rank], &z_status[rank]);
    }


    /*
    std::cout << "rank : " << rank << std::endl;
    for (int i = 0; i < sqrt_p; ++i) {
        for(int j = 0; j< sqrt_p; ++j) {
             std::cout << sub_matrix_X[i][j] << " " << sub_matrix_Y[i][j] << " " << sub_matrix_Z[i][j] << std::endl;
        }
        std::cout << "\n";
    }
    */
    

    int tag = 99;
    int p_src_row = rank / sqrt_p; 
    int p_src_col = rank % sqrt_p;

    int x_send_to_p_dst_row = p_src_row; 
    int x_send_to_p_dst_col = (sqrt_p + p_src_col - p_src_row) % sqrt_p; 
    int x_p_dst_rank = (sqrt_p * x_send_to_p_dst_row ) + x_send_to_p_dst_col;

    int x_recv_from_p_src_row = p_src_row;
    int x_recv_from_p_src_col = (sqrt_p + p_src_col + p_src_row) % sqrt_p; 
    int x_p_src_rank = (sqrt_p * x_recv_from_p_src_row ) + x_recv_from_p_src_col;

    int y_send_to_p_dst_row = (sqrt_p + p_src_row - p_src_col) % sqrt_p; 
    int y_send_to_p_dst_col = p_src_col; 
    int y_p_dst_rank = (sqrt_p * y_send_to_p_dst_row ) + y_send_to_p_dst_col;

    int y_recv_from_p_src_row = (sqrt_p + p_src_row + p_src_col) % sqrt_p;
    int y_recv_from_p_src_col = p_src_col; 
    int y_p_src_rank = (sqrt_p * y_recv_from_p_src_row ) + y_recv_from_p_src_col;

    MPI_Sendrecv_replace(&(sub_matrix_X[0][0]), proc_mat_size * proc_mat_size, MPI_INT,
                         x_p_dst_rank, tag, x_p_src_rank, tag, MPI_COMM_WORLD, &status);

    MPI_Sendrecv_replace(&(sub_matrix_Y[0][0]), proc_mat_size * proc_mat_size, MPI_INT,
                         y_p_dst_rank, tag, y_p_src_rank, tag, MPI_COMM_WORLD, &status);

    for(int l = 0; l < sqrt_p; l++) {

        Matrix_Multiply(sub_matrix_X, sub_matrix_Y, sub_matrix_Z,
                        0, 0, 0, 0, 0, 0, proc_mat_size);

        x_send_to_p_dst_col = (sqrt_p + p_src_col - 1) % sqrt_p; 
        x_p_dst_rank = (sqrt_p * x_send_to_p_dst_row ) + x_send_to_p_dst_col;

        x_recv_from_p_src_col = (sqrt_p + p_src_col + 1) % sqrt_p; 
        x_p_src_rank = (sqrt_p * x_recv_from_p_src_row ) + x_recv_from_p_src_col;

        y_send_to_p_dst_row = (sqrt_p + p_src_row - 1) % sqrt_p; 
        y_send_to_p_dst_col = p_src_col; 
        y_p_dst_rank = (sqrt_p * y_send_to_p_dst_row ) + y_send_to_p_dst_col;

        y_recv_from_p_src_row = (sqrt_p + p_src_row + 1) % sqrt_p;
        y_recv_from_p_src_col = p_src_col; 
        y_p_src_rank = (sqrt_p * y_recv_from_p_src_row ) + y_recv_from_p_src_col;

        MPI_Sendrecv_replace(&(sub_matrix_X[0][0]), proc_mat_size * proc_mat_size, MPI_INT,
                                 x_p_dst_rank, tag, x_p_src_rank, tag, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(&(sub_matrix_Y[0][0]), proc_mat_size * proc_mat_size, MPI_INT,
                                 y_p_dst_rank, tag, y_p_src_rank, tag, MPI_COMM_WORLD, &status);
    }

    int rcv_tag = 232;
    MPI_Status rcv_status[n];
    MPI_Request rcv_request[n];
    /*for (int i = 0; i < sqrt_p; ++i) {
        rcv_tag[i] = 700 + i;
        //std::cout << " receive tag - " << rcv_tag[i] << std::endl;
    }*/
    

    /* 
    std::cout << "before send rank : " << rank << std::endl;
    for (int i = 0; i < proc_mat_size; ++i) {
        for (int j = 0; j < proc_mat_size; ++j) {
            std::cout << sub_matrix_Z[i][j] << "  mmmm ";
        }
        std::cout << "\n";
    }
    */
    
    if (rank != 0) {
	//MPI_Isend(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, rcv_tag, MPI_COMM_WORLD, &rcv_request[rank]);
	MPI_Send(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, rcv_tag, MPI_COMM_WORLD);
    }

    if (rank == 0) {
	update_result(&(sub_matrix_Z[0][0]), 0, Z, n , p);
        for (int i = 1; i < sqrt_p * sqrt_p; ++i) {
            int *temp = (int *)malloc(proc_mat_size*proc_mat_size*sizeof(int));
            //MPI_Irecv(temp, proc_mat_size * proc_mat_size, MPI_INT, i, rcv_tag, MPI_COMM_WORLD, &rcv_request[i]);
            MPI_Recv(temp, proc_mat_size * proc_mat_size, MPI_INT, i, rcv_tag, MPI_COMM_WORLD, &rcv_status[i]);
            //MPI_Wait(&rcv_request[i], &rcv_status[i]);
	    //rcv_request[i] = MPI_REQUEST_NULL;
            /*
            std::cout << "rank : " << i  << " TAG - "<< rcv_tag[i] << std::endl;
            for (int j = 0; j < proc_mat_size * proc_mat_size; ++j) {
                std::cout << temp[j] << " ";
            }
            std::cout << "\n";
            */
            update_result(temp, i, Z, n , p);
            free(temp);
        }
    }


    //MPI_Gatherv(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size,  MPI_INT, mat_Z, send_mat_count, 
    //            send_mat_start, sub_mat_type, 0, MPI_COMM_WORLD);

    /*
    if(rank == 0) {
        std::cout << "\nMASTER OUTPUT:\n";
        printMatrix(Z, n);
    }
    */

    //MPI_TYPE_FREE(&sub_mat_type);
    delete_2d_array(&sub_matrix_X);
    delete_2d_array(&sub_matrix_Y);
    delete_2d_array(&sub_matrix_Z);
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
        //printMatrix(X, n);

        create_2d_array(&Y, n, n);
        fillMatrix(Y, n);

        create_2d_array(&Z, n, n);
        init_sub_matrix(Z, n, n);
    }

    using namespace std::chrono;    
    high_resolution_clock::time_point start_time = high_resolution_clock::now();

    MM_rotate_A_rotate_B(n, processors, myrank);

    high_resolution_clock::time_point end_time = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end_time - start_time);

    if (myrank == 0) {
       delete_2d_array(&X);
       delete_2d_array(&Y);
       delete_2d_array(&Z);
    }

    if(myrank == 0) {
        std::cout << "Exectution Time: " << time_span.count() << " seconds.";
        std::cout << std::endl;
    }
 
    MPI_Finalize();
    return 1;
}
