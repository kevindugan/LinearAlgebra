#include <mpi.h>

class LinearAlgebra {

    public:
        LinearAlgebra(int* argc, char*** argv);
        virtual ~LinearAlgebra();

        int rank() const {return this->comm_rank;}
        int size() const {return this->comm_size;}

    private:
        int comm_rank, comm_size;
};