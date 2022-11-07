
#include "Nucleus.h"
#include "GaussianElimination.h"
#include <random>

#include "Matrix_CycleRowPartition.h"
#include "Vector_CyclePartition.h"

TEST(GaussianElimination, cycle_partition){
    Parallel init;
    const int size = 23;
    Matrix_CycleRowPartition A(size, size, init);
    Vector_CyclePartition x(size, init);
    std::unique_ptr<AbstractVector> b = std::make_unique<Vector_CyclePartition>(size, init);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 10.0);

    std::vector<std::vector<double>> m_vals(size, std::vector<double>(size));
    std::vector<double> v_vals(size);
    for (unsigned int i = 0; i < m_vals.size(); i++){
        for (unsigned int j = 0; j < m_vals.size(); j++)
            m_vals[i][j] = distribution(generator);
        v_vals[i] = distribution(generator);
    }

    A.setValues(m_vals);
    x.setValues(v_vals);

    b = A.mult(x);

    GaussianElimination solver(init);
    std::unique_ptr<AbstractVector> result = solver.solve(A, *b);

    double diff = (x.add(-1.0, *result))->l2norm();
    EXPECT_NEAR(diff, 0.0, 5.0e-11);
}
