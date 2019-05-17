
#include "Nucleus.h"
#include "Matrix.h"

TEST(ParallelMatrix, init_modulo){
    LinearAlgebra init;
    Matrix m(15, 15, init);

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

TEST(ParallelMatrix, init_non_modulo){
  LinearAlgebra init;
  Matrix m(16, 16, init);

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

  Matrix w(128, 128, init);

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

TEST(ParallelMatrix, zero_entries){
  LinearAlgebra init;
  Matrix m(128, 128, init);

  m.setValues(0.12);
  ASSERT_NEAR(m.frobeniusNorm(), 15.36, 1.0e-13);

  m.zeros();
  ASSERT_DOUBLE_EQ(m.frobeniusNorm(), 0.0);
}

TEST(ParallelMatrix, set){
  LinearAlgebra init;
  Matrix m(23, 23, init);

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

TEST(ParallelMatrix, mat_vec_mult){
  LinearAlgebra init;
  Matrix m(6, 6, init);
  Vector v(6, init);

  std::vector<std::vector<double>> m_vals = {{4.3,   0.052, 2.3,   3.1,   0.042, 5.1},
                                             {8.2,   3.6,   4.3,   8.2,   3.6,   4.3},
                                             {0.032, 2.3,   3.1,   0.042, 5.1,   3.1},
                                             {2.3,   3.1,   0.042, 5.1,   6.3,   0.042},
                                             {2.3,   3.1,   0.042, 5.1,   3.1,   4.3},
                                             {0.027, 0.032, 6.3,   8.2,   3.6,   1.4}};

                std::vector<double> v_vals = {5.1,   0.027, 0.027, 0.032, 2.3,   3.1};

  std::vector<double> expected_v = {37.999304, 63.9057, 21.650344, 26.598234, 32.438034, 13.191064};

  Vector expected(6, init);

  m.setValues(m_vals);
  v.setValues(v_vals);
  expected.setValues(expected_v);

  // Perform Mat Vec Mult
  Vector result = m.mult(v);

  ASSERT_EQ(result.size(), expected_v.size());
  for (unsigned int i = 0; i < expected_v.size(); i++)
    ASSERT_DOUBLE_EQ(result.getValue(i), expected_v[i]);

  Vector diff = result.add(-1.0, expected);
  EXPECT_NEAR(diff.length(), 0.0, 1.0E-12);
}
