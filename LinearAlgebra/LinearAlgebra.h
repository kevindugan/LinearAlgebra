#include <mpi.h>

class LinearAlgebra {

    public:
        LinearAlgebra();
        LinearAlgebra(int& argc, char**& argv);
        virtual ~LinearAlgebra();

        int rank() const {return this->comm_rank;}
        int size() const {return this->comm_size;}
        const MPI::Comm& comm() const {return this->communicator;}

    private:
        int comm_rank, comm_size;
        bool linalg_initialized_mpi;
        MPI::Comm& communicator;
};

#define linalg_assert(condition) \
    if(!(condition)){ \
        std::cerr << "Error at " << __FILE__ << ": " << __LINE__ << std::endl; \
        MPI::COMM_WORLD.Abort(1); \
    }

#define info(msg) \
    std::cout << "\e[0;32m[INFO]: \e[0m" << msg << std::endl;
