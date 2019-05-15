
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../Vector.h"

using ::testing::ElementsAreArray;

TEST(ParallelVector, init_modulo){
  LinearAlgebra init;
  Vector v(15, init);

  ASSERT_EQ(v.size(), 15);
  ASSERT_DOUBLE_EQ(v.length(), 0.0);
  std::vector<unsigned int> part = v.getPartitionSize();
  ASSERT_THAT(part, ElementsAreArray({5, 5, 5}));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI::Init(argc, argv);
  int result = RUN_ALL_TESTS();
  MPI::Finalize();
  return result;
}