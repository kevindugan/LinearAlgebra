
#include "Nucleus.h"
#include "Vector.h"

TEST(ParallelVector, init_modulo){
  LinearAlgebra init;
  Vector v(15, init);

  ASSERT_EQ(v.size(), 15);
  ASSERT_DOUBLE_EQ(v.length(), 0.0);
  std::vector<int> part = v.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({5, 5, 5}));

  IndexRange range = v.getGlobalIndexRange();
  std::vector<unsigned int> expected(2);
  if (init.rank() == 0)
    expected = {0, 5};
  else if (init.rank() == 1)
    expected = {5, 10};
  else if (init.rank() == 2)
    expected = {10, 15};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));
}

TEST(ParallelVector, init_non_modulo){
  LinearAlgebra init;
  Vector v(16, init);

  ASSERT_EQ(v.size(), 16);
  ASSERT_DOUBLE_EQ(v.length(), 0.0);
  auto part = v.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({6, 5, 5}));

  IndexRange range = v.getGlobalIndexRange();
  std::vector<unsigned int> expected(2);
  if (init.rank() == 0)
    expected = {0, 6};
  else if (init.rank() == 1)
    expected = {6, 11};
  else if (init.rank() == 2)
    expected = {11, 16};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));

  Vector w(128, init);

  ASSERT_EQ(w.size(), 128);
  ASSERT_DOUBLE_EQ(w.length(), 0.0);
  part = w.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({43, 43, 42}));

  range = w.getGlobalIndexRange();
  if (init.rank() == 0)
    expected = {0, 43};
  else if (init.rank() == 1)
    expected = {43, 86};
  else if (init.rank() == 2)
    expected = {86, 128};

  ASSERT_THAT(expected, ElementsAreArray({range.begin, range.end}));
}

TEST(ParallelVector, zero_entries){
  LinearAlgebra init;
  Vector v(128, init);

  v.setValues(0.12);
  ASSERT_DOUBLE_EQ(v.length(), 1.3576450198781713);
  v.zeros();
  ASSERT_DOUBLE_EQ(v.length(), 0.0);
}

TEST(ParallelVector, set){
  LinearAlgebra init;
  Vector v(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  double expected_norm = 0.0;
  for (auto item : vals)
    expected_norm += item * item;
  expected_norm = sqrt(expected_norm);

  v.setValues(vals);
  ASSERT_DOUBLE_EQ(v.length(), expected_norm);
}

TEST(ParallelVector, getRankWithIndex){
  LinearAlgebra init;
  Vector v(8, init);

  std::vector<unsigned int> ranks = {0, 0, 0, 1, 1, 1, 2, 2};
  for (unsigned int i = 0; i < v.size(); i++)
    ASSERT_EQ(v.findRankWithIndex(i), ranks[i]);
}

TEST(ParallelVector, get){
  LinearAlgebra init;
  Vector v(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  v.setValues(vals);

  for (unsigned int i = 0; i < vals.size(); i++)
    ASSERT_DOUBLE_EQ(v.getValue(i), vals[i]);
}

TEST(ParallelVector, get_global){
  LinearAlgebra init;
  Vector v(15, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2, 9.2, 4.7};
  v.setValues(vals);

  std::vector<double> result(v.size());
  v.getGlobalValues(result);

  ASSERT_EQ(result.size(), vals.size());
  for (unsigned int i = 0; i < result.size(); i++)
    ASSERT_DOUBLE_EQ(result[i], vals[i]);
}

TEST(ParallelVector, copy){
  LinearAlgebra init;
  Vector source(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  source.setValues(vals);

  Vector result1 = source;
  ASSERT_EQ(result1.size(), source.size());
  for (unsigned int i = 0; i < result1.size(); i++)
    EXPECT_DOUBLE_EQ(result1.getValue(i), vals[i]);

  Vector result2(source);
  ASSERT_EQ(result2.size(), source.size());
  for (unsigned int i = 0; i < result2.size(); i++)
    EXPECT_DOUBLE_EQ(result2.getValue(i), vals[i]);

  Vector result3(1, init);
  result3 = source;
  ASSERT_EQ(result3.size(), source.size());
  for (unsigned int i = 0; i < result3.size(); i++)
    EXPECT_DOUBLE_EQ(result3.getValue(i), vals[i]);

}

TEST(ParallelVector, add){
  LinearAlgebra init;
  Vector one(13, init);
  Vector two(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  one.setValues(vals);
  two.setValues(vals);

  Vector result = one.add(two);
  for (unsigned int i = 0; i < vals.size(); i++)
    ASSERT_DOUBLE_EQ(result.getValue(i), 2.0 * vals[i]);

  result = one.add(-1.0, two);
  for (unsigned int i = 0; i < vals.size(); i++)
    ASSERT_DOUBLE_EQ(result.getValue(i), 0.0);

  result = one.add(4.1, two);
  for (unsigned int i = 0; i < vals.size(); i++)
    ASSERT_DOUBLE_EQ(result.getValue(i), vals[i] + 4.1 * vals[i]);
}
