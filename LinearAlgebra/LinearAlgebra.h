#ifndef LINEAR_ALGEBRA_H_A72K
#define LINEAR_ALGEBRA_H_A72K

#include <mpi.h>
#include <iostream>

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
        MPI_Abort(MPI_COMM_WORLD, 1); \
    }

#define info(msg) \
    std::cout << "\e[0;32m[INFO]: \e[0m" << msg << std::endl;

#define Nucleus_ASSERT_EQ(left, right) \
    if(!(left == right)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << " --- " << left << " != " << right << std::endl; \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    }

#endif // LINEAR_ALGEBRA_H_A72K
