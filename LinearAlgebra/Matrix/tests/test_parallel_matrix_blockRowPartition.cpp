
#include "Nucleus.h"
#include "Matrix_BlockRowPartition.h"

TEST(ParallelMatrix_BlockRowPartition, init_modulo){
    LinearAlgebra init;
    Matrix_BlockRowPartition m(15, 15, init);

    ASSERT_EQ(m.nRows(), 15);
    ASSERT_EQ(m.nCols(), 15);
    std::vector<unsigned int> part = m.getPartitionSize();
    ASSERT_THAT(part, ElementsAreArray({3, 3, 3, 3, 3}));

    IndexRange range = m.getGlobalRowIndexRange();
    std::vector<unsigned int> expected(2);
    if (init.rank() == 0)
        expected = {0, 3};
    else if (init.rank() == 1)
        expected = {3, 6};
    else if (init.rank() == 2)
        expected = {6, 9};
    else if (init.rank() == 3)
        expected = {9, 12};
    else if (init.rank() == 4)
        expected = {12, 15};

    ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));
}

TEST(ParallelMatrix_BlockRowPartition, init_non_modulo){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(16, 16, init);

  ASSERT_EQ(m.nRows(), 16);
  ASSERT_EQ(m.nCols(), 16);
  auto part = m.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({4, 3, 3, 3, 3}));

  IndexRange range = m.getGlobalRowIndexRange();
  std::vector<unsigned int> expected(2);
  if (init.rank() == 0)
    expected = {0, 4};
  else if (init.rank() == 1)
    expected = {4, 7};
  else if (init.rank() == 2)
    expected = {7, 10};
  else if (init.rank() == 3)
    expected = {10, 13};
  else if (init.rank() == 4)
    expected = {13, 16};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));

  Matrix_BlockRowPartition w(128, 128, init);

  ASSERT_EQ(w.nRows(), 128);
  ASSERT_EQ(w.nCols(), 128);
  part = w.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({26, 26, 26, 25, 25}));

  range = w.getGlobalRowIndexRange();
  if (init.rank() == 0)
    expected = {0, 26};
  else if (init.rank() == 1)
    expected = {26, 52};
  else if (init.rank() == 2)
    expected = {52, 78};
  else if (init.rank() == 3)
    expected = {78, 103};
  else if (init.rank() == 4)
    expected = {103, 128};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));
}

TEST(ParallelMatrix_BlockRowPartition, zero_entries){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(128, 128, init);

  m.setValues(0.12);
  ASSERT_NEAR(m.frobeniusNorm(), 15.36, 1.0e-13);

  m.zeros();
  ASSERT_DOUBLE_EQ(m.frobeniusNorm(), 0.0);
}

TEST(ParallelMatrix_BlockRowPartition, set){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(23, 23, init);

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

TEST(ParallelMatrix_BlockRowPartition, mat_vec_mult){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(6, 6, init);
  Vector_BlockPartition v(6, init);

  std::vector<std::vector<double>> m_vals = {{4.3,   0.052, 2.3,   3.1,   0.042, 5.1},
                                             {8.2,   3.6,   4.3,   8.2,   3.6,   4.3},
                                             {0.032, 2.3,   3.1,   0.042, 5.1,   3.1},
                                             {2.3,   3.1,   0.042, 5.1,   6.3,   0.042},
                                             {2.3,   3.1,   0.042, 5.1,   3.1,   4.3},
                                             {0.027, 0.032, 6.3,   8.2,   3.6,   1.4}};

                std::vector<double> v_vals = {5.1,   0.027, 0.027, 0.032, 2.3,   3.1};

  std::vector<double> expected_v = {37.999304, 63.9057, 21.650344, 26.598234, 32.438034, 13.191064};

  Vector_BlockPartition expected(6, init);

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

TEST(ParallelMatrix_BlockRowPartition, getRankWithIndex){
  LinearAlgebra init;
  Matrix_BlockRowPartition v(8, 8, init);

  std::vector<unsigned int> ranks = {0, 0, 1, 1, 2, 2, 3, 4};
  for (unsigned int i = 0; i < v.nRows(); i++)
    ASSERT_EQ(v.findRankWithIndex(i), ranks[i]);
}

TEST(ParallelMatrix_BlockRowPartition, get){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(23, 23, init);

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

TEST(ParallelMatrix_BlockRowPartition, copy){
    LinearAlgebra init;
    Matrix_BlockRowPartition source(23, 23, init);

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

    Matrix_BlockRowPartition result1 = source;
    ASSERT_EQ(result1.nRows(), vals.size());
    for (unsigned int row = 0; row < result1.nRows(); row++){
        ASSERT_EQ(result1.nCols(), vals[row].size());
        for (unsigned int col = 0; col < result1.nCols(); col++)
            EXPECT_DOUBLE_EQ(result1.getValue(row,col), vals[row][col]);
    }

    Matrix_BlockRowPartition result2(source);
    ASSERT_EQ(result2.nRows(), vals.size());
    for (unsigned int row = 0; row < result2.nRows(); row++){
        ASSERT_EQ(result2.nCols(), vals[row].size());
        for (unsigned int col = 0; col < result2.nCols(); col++)
            EXPECT_DOUBLE_EQ(result2.getValue(row,col), vals[row][col]);
    }

    Matrix_BlockRowPartition result3(1, 1, init);
    result3 = source;
    ASSERT_EQ(result3.nRows(), vals.size());
    for (unsigned int row = 0; row < result3.nRows(); row++){
        ASSERT_EQ(result3.nCols(), vals[row].size());
        for (unsigned int col = 0; col < result3.nCols(); col++)
            EXPECT_DOUBLE_EQ(result3.getValue(row,col), vals[row][col]);
    }
}

TEST(ParallelMatrix_BlockRowPartition, clone){
    LinearAlgebra init;
    Matrix_BlockRowPartition source(23, 23, init);

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

TEST(ParallelMatrix_BlockRowPartition, get_row){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(23, 23, init);

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
  std::vector<unsigned int> offset(init.size());
  for (unsigned int i = 0, sum = 0; i < part.size(); i++, sum+=part[i-1])
    offset[i] = sum;

  for (unsigned int i = 0; i < part[init.rank()]; i++){
    std::vector<double> row_values = m.getLocalRowValues(i);
    for (unsigned int j = 0; j < m.nCols(); j++)
      EXPECT_DOUBLE_EQ(row_values[j], vals[offset[init.rank()] + i][j]);
  }
}

TEST(ParallelMatrix_BlockRowPartition, set_row){
  LinearAlgebra init;
  Matrix_BlockRowPartition m(23, 23, init);

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
  std::vector<unsigned int> offset(init.size());
  for (unsigned int i = 0, sum = 0; i < part.size(); i++, sum+=part[i-1])
    offset[i] = sum;

  m.zeros();

  for (unsigned int row = 0; row < part[init.rank()]; row++)
    m.setLocalRowValues(row, vals[offset[init.rank()] + row]);

  // Check values
  for (unsigned int row = offset[init.rank()]; row < part[init.rank()]; row++){
      for (unsigned int col = 0; col < m.nCols(); col++)
          EXPECT_DOUBLE_EQ(m.getValue(row,col), vals[row][col]);
  }
}