#include "Nucleus.h"
#include "Vector_CyclePartition.h"

TEST(ParallelVector_CyclePartition, inti_modulo){
    LinearAlgebra init;
    Vector_CyclePartition v(16, init);

    ASSERT_EQ(v.size(), 16);
    ASSERT_DOUBLE_EQ(v.l2norm(), 0.0);
    std::vector<int> part = v.getPartitionSize();
    EXPECT_THAT(part, ElementsAreArray({4, 4, 4, 4}));

    IndexRange range = v.getGlobalIndexRange();
    std::vector<unsigned int> expected(3);
    if (init.rank() == 0)
        expected = {0, 16, 4};
    else if (init.rank() == 1)
        expected = {1, 16, 4};
    else if (init.rank() == 2)
        expected = {2, 16, 4};
    else if (init.rank() == 3)
        expected = {3, 16, 4};

    EXPECT_THAT(expected, ElementsAreArray({range.begin, range.end, range.skip}));
}

TEST(ParallelVector_CyclePartition, init_non_modulo){
  LinearAlgebra init;
  Vector_CyclePartition v(17, init);

  ASSERT_EQ(v.size(), 17);
  EXPECT_DOUBLE_EQ(v.l2norm(), 0.0);
  auto part = v.getPartitionSize();
  EXPECT_THAT(part, ElementsAreArray({5, 4, 4, 4}));

  IndexRange range = v.getGlobalIndexRange();
  std::vector<unsigned int> expected(3);
  if (init.rank() == 0)
      expected = {0, 17, 4};
  else if (init.rank() == 1)
      expected = {1, 17, 4};
  else if (init.rank() == 2)
      expected = {2, 17, 4};
  else if (init.rank() == 3)
      expected = {3, 17, 4};

  EXPECT_THAT(expected, ElementsAreArray({range.begin, range.end, range.skip}));

  Vector_CyclePartition w(131, init);

  ASSERT_EQ(w.size(), 131);
  EXPECT_DOUBLE_EQ(w.l2norm(), 0.0);
  part = w.getPartitionSize();
  EXPECT_THAT(part, ElementsAreArray({33, 33, 33, 32}));

  range = w.getGlobalIndexRange();
  if (init.rank() == 0)
      expected = {0, 131, 4};
  else if (init.rank() == 1)
      expected = {1, 131, 4};
  else if (init.rank() == 2)
      expected = {2, 131, 4};
  else if (init.rank() == 3)
      expected = {3, 131, 4};

  EXPECT_THAT(expected, ElementsAreArray({range.begin, range.end, range.skip}));
}

TEST(ParallelVector_CyclePartition, zero_entries){
  LinearAlgebra init;
  Vector_CyclePartition v(131, init);

  v.setValues(0.12);
  EXPECT_DOUBLE_EQ(v.l2norm(), sqrt( 131.0 * 0.12 * 0.12));
  v.zeros();
  EXPECT_DOUBLE_EQ(v.l2norm(), 0.0);
}

TEST(ParallelVector_CyclePartition, set){
  LinearAlgebra init;
  Vector_CyclePartition v(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  double expected_norm = 0.0;
  for (auto item : vals)
    expected_norm += item * item;
  expected_norm = sqrt(expected_norm);

  v.setValues(vals);
  EXPECT_DOUBLE_EQ(v.l2norm(), expected_norm);
}

TEST(ParallelVector_CyclePartition, getRankWithIndex){
  LinearAlgebra init;
  Vector_CyclePartition v(10, init);

  std::vector<unsigned int> ranks = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1};
  for (unsigned int i = 0; i < v.size(); i++)
    EXPECT_EQ(v.findRankWithIndex(i), ranks[i]);
}

TEST(ParallelVector_CyclePartition, get){
  LinearAlgebra init;
  Vector_CyclePartition v(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  v.setValues(vals);

  for (unsigned int i = 0; i < vals.size(); i++)
    EXPECT_DOUBLE_EQ(v.getValue(i), vals[i]);
}

TEST(ParallelVector_CyclePartition, get_global){
  LinearAlgebra init;
  Vector_CyclePartition v(15, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2, 9.2, 4.7};
  v.setValues(vals);

  std::vector<double> result(v.size());
  v.getGlobalValues(result);

  ASSERT_EQ(result.size(), vals.size());
  for (unsigned int i = 0; i < result.size(); i++)
    EXPECT_DOUBLE_EQ(result[i], vals[i]);
}

TEST(ParallelVector_CyclePartition, copy){
  LinearAlgebra init;
  Vector_CyclePartition source(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  source.setValues(vals);

  Vector_CyclePartition result1 = source;
  ASSERT_EQ(result1.size(), source.size());
  for (unsigned int i = 0; i < result1.size(); i++)
    EXPECT_DOUBLE_EQ(result1.getValue(i), vals[i]);

  Vector_CyclePartition result2(source);
  ASSERT_EQ(result2.size(), source.size());
  for (unsigned int i = 0; i < result2.size(); i++)
    EXPECT_DOUBLE_EQ(result2.getValue(i), vals[i]);

  Vector_CyclePartition result3(1, init);
  result3 = source;
  ASSERT_EQ(result3.size(), source.size());
  for (unsigned int i = 0; i < result3.size(); i++)
    EXPECT_DOUBLE_EQ(result3.getValue(i), vals[i]);
}

TEST(ParallelVector_CyclePartition, clone){
  LinearAlgebra init;
  Vector_CyclePartition source(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  source.setValues(vals);

  std::unique_ptr<AbstractVector> result1 = source.clone();
  ASSERT_EQ(result1->size(), source.size());
  for (unsigned int i = 0; i < result1->size(); i++)
    EXPECT_DOUBLE_EQ(result1->getValue(i), vals[i]);
}

TEST(ParallelVector_CyclePartition, add){
  LinearAlgebra init;
  Vector_CyclePartition one(13, init);
  Vector_CyclePartition two(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  one.setValues(vals);
  two.setValues(vals);

  std::unique_ptr<AbstractVector> result = one.add(two);
  for (unsigned int i = 0; i < vals.size(); i++)
    EXPECT_DOUBLE_EQ(result->getValue(i), 2.0 * vals[i]);

  result = one.add(-1.0, two);
  for (unsigned int i = 0; i < vals.size(); i++)
    EXPECT_DOUBLE_EQ(result->getValue(i), 0.0);

  result = one.add(4.1, two);
  for (unsigned int i = 0; i < vals.size(); i++)
    EXPECT_DOUBLE_EQ(result->getValue(i), vals[i] + 4.1 * vals[i]);
}

TEST(ParallelVector_CyclePartition, print){
  LinearAlgebra init;
  Vector_CyclePartition v(13, init);

  std::vector<double> vals = {1.1, 2.3, 3.1, 4.2, 5.1, 2.7, 3.2, 6.3, 8.2, 3.6, 4.3, 5.2, 7.2};
  v.setValues(vals);

  std::stringstream expected;
  for (unsigned int i = 0; i < vals.size(); i++)
    expected << std::setw(9) << std::setprecision(2) << std::scientific << vals[i] << std::endl;

  std::stringstream result;
  v.print(result);

  if (init.rank() == 0){
    ASSERT_EQ(expected.str().size(), result.str().size());
    EXPECT_STREQ(result.str().c_str(), expected.str().c_str());
  }

}
