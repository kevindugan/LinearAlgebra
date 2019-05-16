
#include "Nucleus.h"
#include "Matrix.h"

TEST(ParallelMatrix, init_modulo){
    LinearAlgebra init;
    Matrix m(15, 15, init);

    ASSERT_EQ(m.nRows(), 15);
    ASSERT_EQ(m.nCols(), 15);
    std::vector<unsigned int> part = m.getPartitionSize();
    ASSERT_THAT(part, ElementsAreArray({5, 5, 5}));

    IndexRange range = m.getGlobalRowIndexRange();
    std::vector<unsigned int> expected(2);
    if (init.rank() == 0)
        expected = {0, 5};
    else if (init.rank() == 1)
        expected = {5, 10};
    else if (init.rank() == 2)
        expected = {10, 15};

    ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));
}

TEST(ParallelMatrix, init_non_modulo){
  LinearAlgebra init;
  Matrix m(16, 16, init);

  ASSERT_EQ(m.nRows(), 16);
  ASSERT_EQ(m.nCols(), 16);
  auto part = m.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({6, 5, 5}));

  IndexRange range = m.getGlobalRowIndexRange();
  std::vector<unsigned int> expected(2);
  if (init.rank() == 0)
    expected = {0, 6};
  else if (init.rank() == 1)
    expected = {6, 11};
  else if (init.rank() == 2)
    expected = {11, 16};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));

  Matrix w(128, 128, init);

  ASSERT_EQ(w.nRows(), 128);
  ASSERT_EQ(w.nCols(), 128);
  part = w.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({43, 43, 42}));

  range = w.getGlobalRowIndexRange();
  if (init.rank() == 0)
    expected = {0, 43};
  else if (init.rank() == 1)
    expected = {43, 86};
  else if (init.rank() == 2)
    expected = {86, 128};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));
}

TEST(ParallelMatrix, zero_entries){
  LinearAlgebra init;
  Matrix m(128, 128, init);

//   m.setValues(0.12);
//   ASSERT_NE(m.length(), 0.0);
//   m.zeros();
//   ASSERT_DOUBLE_EQ(m.length(), 0.0);
}