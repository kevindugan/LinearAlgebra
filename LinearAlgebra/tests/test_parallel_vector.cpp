
#include "gtest/gtest.h"
#include "../Vector.h"

TEST(ParallelVector, norm){
    LinearAlgebra init;
    Vector v(15, init);

    ASSERT_EQ(v.size(), 15);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  int result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}