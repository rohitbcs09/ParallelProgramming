#include <iostream>
#include <stdio.h>
#include <mpi.h>

int main( int argc, char *argv[] ){
    std::cout << "Argc = " << argc<< " " << argv[1] << "\n";
    int myrank, v = 121;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (myrank == 0) {
        MPI_Send(&v, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("Process %d sent %d \n", myrank, v );
    } 
    else if (myrank == 1) {
        MPI_Recv(&v, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD , &status );
        printf("Process %d received %d \n", myrank, v );
    }
    MPI_Finalize();
}
