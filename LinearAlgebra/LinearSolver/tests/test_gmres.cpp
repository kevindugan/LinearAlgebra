
#include "Nucleus.h"
#include "GMRes.h"
#include <random>

#include "SparseMatrix_BlockRowPartition.h"
#include "Vector_BlockPartition.h"

TEST(GMRes, sparseBlock){
    GTEST_SKIP() << " GMRes not yet implemented";

    Parallel init;
    const int size = 23;
    SparseMatrix_BlockRowPartition A(size, size, init);
    Vector_BlockPartition x(size, init);

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

    GMRes solver(init);
    std::unique_ptr<AbstractVector> result = solver.solve(A, *b);

    double diff = (x.add(-1.0, *result))->l2norm();
    EXPECT_NEAR(diff, 0.0, 5.0e-11);
}
