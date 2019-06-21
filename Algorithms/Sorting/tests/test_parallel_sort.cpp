#include <parallel_sorting>
#include <numeric>
#include "Nucleus.h"
#include "Parallel.h"

TEST(Parallel_Sort, no_overlap){
    Parallel init;
    std::vector<unsigned int> vec(10);
    std::iota(vec.begin(), vec.end(), 0);
}