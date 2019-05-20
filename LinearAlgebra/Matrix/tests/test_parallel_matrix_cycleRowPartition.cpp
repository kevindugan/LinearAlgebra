
#include "Nucleus.h"
#include "Matrix_CycleRowPartition.h"

TEST(ParallelMatrix_CycleRowPartition, init_modulo){
    LinearAlgebra init;
    Matrix_CycleRowPartition m(16, 16, init);

    ASSERT_EQ(m.nRows(), 16);
    ASSERT_EQ(m.nCols(), 16);
    std::vector<unsigned int> part = m.getPartitionSize();
    ASSERT_THAT(part, ElementsAreArray({4, 4, 4, 4}));

    IndexRange range = m.getGlobalRowIndexRange();
    std::vector<unsigned int> expected(3);
    if (init.rank() == 0)
        expected = {0, 16, 4};
    else if (init.rank() == 1)
        expected = {1, 16, 4};
    else if (init.rank() == 2)
        expected = {2, 16, 4};
    else if (init.rank() == 3)
        expected = {3, 16, 4};

    ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end, range.skip}));
}

TEST(ParallelMatrix_CycleRowPartition, init_non_modulo){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(17, 17, init);

  ASSERT_EQ(m.nRows(), 17);
  ASSERT_EQ(m.nCols(), 17);
  auto part = m.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({5, 4, 4, 4}));

  IndexRange range = m.getGlobalRowIndexRange();
  std::vector<unsigned int> expected(3);
    if (init.rank() == 0)
        expected = {0, 17, 4};
    else if (init.rank() == 1)
        expected = {1, 17, 4};
    else if (init.rank() == 2)
        expected = {2, 17, 4};
    else if (init.rank() == 3)
        expected = {3, 17, 4};

    ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end, range.skip}));

  Matrix_CycleRowPartition w(131, 131, init);

  ASSERT_EQ(w.nRows(), 131);
  ASSERT_EQ(w.nCols(), 131);
  part = w.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({33, 33, 33, 32}));

  range = w.getGlobalRowIndexRange();
  if (init.rank() == 0)
    expected = {0, 131, 4};
  else if (init.rank() == 1)
    expected = {1, 131, 4};
  else if (init.rank() == 2)
    expected = {2, 131, 4};
  else if (init.rank() == 3)
    expected = {3, 131, 4};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end, range.skip}));
}

TEST(ParallelMatrix_CycleRowPartition, zero_entries){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(131, 131, init);

  m.setValues(0.12);
  ASSERT_NEAR(m.frobeniusNorm(), sqrt(0.12 * 0.12 * 131 * 131), 1.0e-12);

  m.zeros();
  ASSERT_DOUBLE_EQ(m.frobeniusNorm(), 0.0);
}

TEST(ParallelMatrix_CycleRowPartition, set){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(23, 23, init);

  std::vector<std::vector<double>> vals = {{1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6}};

  double expected_norm = 0.0;
  for (const auto &row : vals)
    for (const auto &item : row)
      expected_norm += item * item;
  expected_norm = sqrt(expected_norm);

  m.setValues(vals);
  ASSERT_DOUBLE_EQ(m.frobeniusNorm(), expected_norm);
}

TEST(ParallelMatrix_CycleRowPartition, mat_vec_mult){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(6, 6, init);
  Vector_CyclePartition v(6, init);

  std::vector<std::vector<double>> m_vals = {{4.3,   0.052, 2.3,   3.1,   0.042, 5.1},
                                             {8.2,   3.6,   4.3,   8.2,   3.6,   4.3},
                                             {0.032, 2.3,   3.1,   0.042, 5.1,   3.1},
                                             {2.3,   3.1,   0.042, 5.1,   6.3,   0.042},
                                             {2.3,   3.1,   0.042, 5.1,   3.1,   4.3},
                                             {0.027, 0.032, 6.3,   8.2,   3.6,   1.4}};

                std::vector<double> v_vals = {5.1,   0.027, 0.027, 0.032, 2.3,   3.1};

  std::vector<double> expected_v = {37.999304, 63.9057, 21.650344, 26.598234, 32.438034, 13.191064};

  Vector_CyclePartition expected(6, init);

  m.setValues(m_vals);
  v.setValues(v_vals);
  expected.setValues(expected_v);

  // Perform Mat Vec Mult
  std::unique_ptr<AbstractVector> result = m.mult(v);

  ASSERT_EQ(result->size(), expected_v.size());
  for (unsigned int i = 0; i < expected_v.size(); i++)
    ASSERT_DOUBLE_EQ(result->getValue(i), expected_v[i]);

  std::unique_ptr<AbstractVector> diff = result->add(-1.0, expected);
  EXPECT_NEAR(diff->l2norm(), 0.0, 1.0E-12);
}

TEST(ParallelMatrix_CycleRowPartition, print){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(6, 6, init);

  std::vector<std::vector<double>> m_vals = {{4.3,   0.052, 2.3,   3.1,   0.042, 5.1},
                                             {8.2,   3.6,   4.3,   8.2,   3.6,   4.3},
                                             {0.032, 2.3,   3.1,   0.042, 5.1,   3.1},
                                             {2.3,   3.1,   0.042, 5.1,   6.3,   0.042},
                                             {2.3,   3.1,   0.042, 5.1,   3.1,   4.3},
                                             {0.027, 0.032, 6.3,   8.2,   3.6,   1.4}};

    std::stringstream expected;
    for (const auto &line : m_vals){
        for (const auto &item : line)
            expected << std::setw(6) << std::setprecision(4) << item << ", ";
        expected << "\n";
    }

    m.setValues(m_vals);
    std::stringstream result;
    m.print(result);

    if (init.rank() == 0){
        ASSERT_EQ(expected.str().size(), result.str().size());
        EXPECT_STREQ(expected.str().c_str(), result.str().c_str());
    }
}

TEST(ParallelMatrix_CycleRowPartition, get){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(23, 23, init);

  std::vector<std::vector<double>> vals = {{1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6}};

  m.setValues(vals);

  for (unsigned int i = 0; i < m.nRows(); i++)
    for (unsigned int j = 0; j < m.nCols(); j++)
        EXPECT_DOUBLE_EQ(m.getValue(i,j), vals[i][j]);
}

TEST(ParallelMatrix_CycleRowPartition, getRankWithIndex){
  LinearAlgebra init;
  Matrix_CycleRowPartition v(10, 10, init);

  std::vector<unsigned int> ranks = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1};
  for (unsigned int i = 0; i < v.nRows(); i++)
    EXPECT_EQ(v.findRankWithIndex(i), ranks[i]);
}

TEST(ParallelMatrix_CycleRowPartition, copy){
    LinearAlgebra init;
    Matrix_CycleRowPartition source(23, 23, init);

    std::vector<std::vector<double>> vals = {{1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6}};
    source.setValues(vals);

    Matrix_CycleRowPartition result1 = source;
    ASSERT_EQ(result1.nRows(), vals.size());
    for (unsigned int row = 0; row < result1.nRows(); row++){
        ASSERT_EQ(result1.nCols(), vals[row].size());
        for (unsigned int col = 0; col < result1.nCols(); col++)
            EXPECT_DOUBLE_EQ(result1.getValue(row,col), vals[row][col]);
    }

    Matrix_CycleRowPartition result2(source);
    ASSERT_EQ(result2.nRows(), vals.size());
    for (unsigned int row = 0; row < result2.nRows(); row++){
        ASSERT_EQ(result2.nCols(), vals[row].size());
        for (unsigned int col = 0; col < result2.nCols(); col++)
            EXPECT_DOUBLE_EQ(result2.getValue(row,col), vals[row][col]);
    }

    Matrix_CycleRowPartition result3(1, 1, init);
    result3 = source;
    ASSERT_EQ(result3.nRows(), vals.size());
    for (unsigned int row = 0; row < result3.nRows(); row++){
        ASSERT_EQ(result3.nCols(), vals[row].size());
        for (unsigned int col = 0; col < result3.nCols(); col++)
            EXPECT_DOUBLE_EQ(result3.getValue(row,col), vals[row][col]);
    }
}

TEST(ParallelMatrix_CycleRowPartition, clone){
    LinearAlgebra init;
    Matrix_CycleRowPartition source(23, 23, init);

    std::vector<std::vector<double>> vals = {{1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                             {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                             {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                             {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6}};
    source.setValues(vals);

    std::unique_ptr<AbstractMatrix> result1 = source.clone();
    ASSERT_EQ(result1->nRows(), vals.size());
    for (unsigned int row = 0; row < result1->nRows(); row++){
        ASSERT_EQ(result1->nCols(), vals[row].size());
        for (unsigned int col = 0; col < result1->nCols(); col++)
            EXPECT_DOUBLE_EQ(result1->getValue(row,col), vals[row][col]);
    }
}

TEST(ParallelMatrix_CycleRowPartition, get_row){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(23, 23, init);

  std::vector<std::vector<double>> vals = {{1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6}};

  m.setValues(vals);

  // Global Row values
  for (unsigned int i = 0; i < m.nRows(); i++){
    std::vector<double> row_values = m.getRowValues(i);
    for (unsigned int j = 0; j < m.nCols(); j++)
      EXPECT_DOUBLE_EQ(row_values[j], vals[i][j]);
  }

  // Local Row values
  std::vector<unsigned int> part = m.getPartitionSize();
  for (unsigned int i = 0, k=init.rank(); i < part[init.rank()]; i++, k+=init.size()){
    std::vector<double> row_values = m.getLocalRowValues(i);
    for (unsigned int j = 0; j < m.nCols(); j++)
      EXPECT_DOUBLE_EQ(row_values[j], vals[k][j]);
  }
}

TEST(ParallelMatrix_CycleRowPartition, set_row){
  LinearAlgebra init;
  Matrix_CycleRowPartition m(23, 23, init);

  std::vector<std::vector<double>> vals = {{1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 6.3, 0.042, 0.032, 6.3, 8.2, 3.6, 4.3, 8.2, 3.6, 4.3, 0.052, 5.1, 0.027, 0.032, 6.3},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {1.1, 2.3, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 0.052, 2.3, 3.1, 2.3, 3.1, 0.042, 5.1, 0.027, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6},
                                           {3.1, 4.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.6, 3.1, 4.3, 0.052, 7.2, 3.1, 0.042, 5.1, 0.052, 4.3, 0.052, 7.2, 3.1, 0.042, 2.3, 3.1},
                                           {4.3, 0.027, 0.032, 6.3, 8.2, 0.052, 2.3, 0.027, 0.032, 2.3, 3.1, 0.042, 5.1, 3.1, 4.3, 0.052, 2.3, 3.1, 0.042, 0.042, 5.1, 0.027, 0.027},
                                           {5.1, 0.027, 0.032, 6.3, 8.2, 3.6, 4.3, 0.052, 2.3, 3.1, 0.042, 5.1, 0.027, 0.027, 0.032, 6.3, 8.2, 3.6, 0.027, 0.032, 6.3, 8.2, 3.6}};

  EXPECT_DOUBLE_EQ(m.frobeniusNorm(), 0.0);

  // Global Row values
  for (unsigned int i = 0; i < m.nRows(); i++)
    m.setRowValues(i, vals[i]);

  // Check values
  for (unsigned int row = 0; row < m.nRows(); row++){
      ASSERT_EQ(m.nCols(), vals[row].size());
      for (unsigned int col = 0; col < m.nCols(); col++)
          EXPECT_DOUBLE_EQ(m.getValue(row,col), vals[row][col]);
  }

  // Local Row values
  std::vector<unsigned int> part = m.getPartitionSize();
  m.zeros();

  for (unsigned int row = 0, global_row=init.rank(); row < part[init.rank()]; row++, global_row+=init.size())
    m.setLocalRowValues(row, vals[global_row]);

  // Check values
  for (unsigned int row = init.rank(); row < m.nRows(); row+=init.size()){
      for (unsigned int col = 0; col < m.nCols(); col++)
          EXPECT_DOUBLE_EQ(m.getValue(row,col), vals[row][col]);
  }
}