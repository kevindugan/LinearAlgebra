#include "LinearAlgebra.h"

LinearAlgebra::LinearAlgebra(int* argc, char*** argv) {

    MPI_Init(argc, argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &this->comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &this->comm_size);
}

LinearAlgebra::~LinearAlgebra(){
    MPI_Finalize();
}
