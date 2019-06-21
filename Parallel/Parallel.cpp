#include "Parallel.h"

Parallel::Parallel(){
    int initialized;
    MPI_Initialized(&initialized);
    if(initialized){
        MPI_Comm_rank(MPI_COMM_WORLD, &this->comm_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &this->comm_size);
        this->linalg_initialized_mpi = false;
    }
}

Parallel::Parallel(int* argc, char*** argv){
    MPI_Init(argc, argv);
    this->linalg_initialized_mpi = true;
    MPI_Comm_rank(MPI_COMM_WORLD, &this->comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &this->comm_size);
}

Parallel::~Parallel(){
    if (this->linalg_initialized_mpi)
        MPI_Finalize();
}
