
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
