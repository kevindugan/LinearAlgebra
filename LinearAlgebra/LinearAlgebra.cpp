#include "LinearAlgebra.h"

LinearAlgebra::LinearAlgebra() : communicator(MPI::COMM_WORLD){
    if(MPI::Is_initialized()){
        this->comm_rank = this->communicator.Get_rank();
        this->comm_size = this->communicator.Get_size();
        this->linalg_initialized_mpi = false;
    }
}

LinearAlgebra::LinearAlgebra(int& argc, char**& argv) : communicator(MPI::COMM_WORLD) {
    MPI::Init(argc, argv);
    this->linalg_initialized_mpi = true;
    this->comm_rank = this->communicator.Get_rank();
    this->comm_size = this->communicator.Get_size();
}

LinearAlgebra::~LinearAlgebra(){
    if (this->linalg_initialized_mpi)
        MPI::Finalize();
}
