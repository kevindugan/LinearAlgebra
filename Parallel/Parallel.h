#ifndef PARALLEL_H_A72K
#define PARALLEL_H_A72K

#include <mpi.h>
#include <iostream>

class Parallel {

    public:
        Parallel();
        Parallel(int* argc, char*** argv);
        virtual ~Parallel();

        int rank() const {return this->comm_rank;}
        int size() const {return this->comm_size;}

    private:
        int comm_rank, comm_size;
        bool linalg_initialized_mpi;
};

#define Nucleus_ASSERT(condition) \
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

#define Nucleus_ASSERT_LT(left, right) \
    if(!(left < right)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << " --- " << left << " >= " << right << std::endl; \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    }

#define Nucleus_ASSERT_GT(left, right) \
    if(!(left > right)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << " --- " << left << " <= " << right << std::endl; \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    }

#define Nucleus_ASSERT_LE(left, right) \
    if(!(left <= right)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << " --- " << left << " > " << right << std::endl; \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    }

#define Nucleus_ASSERT_GE(left, right) \
    if(!(left >= right)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << " --- " << left << " < " << right << std::endl; \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    }

#endif // PARALLEL_H_A72K
