#include "GaussianElimination.h"
#include <chrono>
#include <random>
#include <iomanip>

#include "Matrix_CycleRowPartition.h"
#include "Vector_CyclePartition.h"

using timer = std::chrono::high_resolution_clock;


int main (int argc, char* argv[]){
    LinearAlgebra init(&argc, &argv);

    // Setup calculation
    timer::time_point start = timer::now();

    // Get system size
    int size = 1e1;
    if (init.rank() == 0){
        for (unsigned int i = 0; i < argc; i++){
            if (std::string(argv[i]) == "--size"){
                size = int(std::atof(argv[i+1]));
            }
        }
        std::cout << "System Size: " << size << std::endl;
    }
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);



    Matrix_CycleRowPartition A(size, size, init);
    Vector_CyclePartition x(size, init);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-10.0, 10.0);

    std::vector<std::vector<double>> m_vals(size, std::vector<double>(size));
    std::vector<double> v_vals(size);
    for (unsigned int i = 0; i < m_vals.size(); i++){
        for (unsigned int j = 0; j < m_vals.size(); j++)
            m_vals[i][j] = distribution(generator);
        v_vals[i] = distribution(generator);
    }

    A.setValues(m_vals);
    x.setValues(v_vals);
    std::unique_ptr<AbstractVector> b = A.mult(x);

    GaussianElimination solver(init);

    double secs = std::chrono::duration_cast<std::chrono::duration<double>>(timer::now()-start).count();
    if (init.rank() == 0){
        std::cout << "[Rank " << std::setw(3) << init.rank() << "] Elapsed time [ init]: " << secs << " s" << std::endl;
        for (unsigned int proc = 1; proc < init.size(); proc++){
            MPI_Recv(&secs, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "[Rank " << std::setw(3) << proc << "] Elapsed time [ init]: " << secs << " s" << std::endl;
        }
    } else {
        MPI_Send(&secs, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    // Perform Calculation
    MPI_Barrier(MPI_COMM_WORLD);
    start = timer::now();
    std::unique_ptr<AbstractVector> result = solver.solve(A, *b);
    secs = std::chrono::duration_cast<std::chrono::duration<double>>(timer::now()-start).count();
    if (init.rank() == 0){
        std::cout << "[Rank " << std::setw(3) << init.rank() << "] Elapsed time [solve]: " << secs << " s" << std::endl;
        for (unsigned int proc = 1; proc < init.size(); proc++){
            MPI_Recv(&secs, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "[Rank " << std::setw(3) << proc << "] Elapsed time [solve]: " << secs << " s" << std::endl;            
        }
    } else {
        MPI_Send(&secs, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    // Check Correctness of solution
    double diff = (x.add(-1.0, *result))->l2norm();
    if (init.rank() == 0){
        std::cout << "Solution Error: " << diff << std::endl;
    }

    return 0;
}