#include <mpi.h>

class LinearAlgebra {

    public:
        LinearAlgebra();
        LinearAlgebra(int* argc, char*** argv);
        virtual ~LinearAlgebra();

        int rank() const {return this->comm_rank;}
        int size() const {return this->comm_size;}

    private:
        int comm_rank, comm_size;
        bool linalg_initialized_mpi;
};

#define linalg_assert(condition) \
    if(!(condition)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << std::endl; \
        MPI_Abort(MPI_COMM_WORLD, 0); \
    }

#define info(msg) \
    std::cout << "\e[0;32m[INFO]: \e[0m" << msg << std::endl;
