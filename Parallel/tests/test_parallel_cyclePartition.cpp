#include "Nucleus.h"
#include "CyclePartition1D.h"

TEST(CyclePartition, init_modulo){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    CyclePartition1D part(15, init);
    ASSERT_EQ(15, part.getGlobalSize());
    std::vector<unsigned int> result = part.getGlobalIndexRange();

    std::vector<unsigned int> expected(3);
    if (init.rank() == 0)
        expected = {0, 15, 5};
    else if (init.rank() == 1)
        expected = {1, 15, 5};
    else if (init.rank() == 2)
        expected = {2, 15, 5};
    else if (init.rank() == 3)
        expected = {3, 15, 5};
    else if (init.rank() == 4)
        expected = {4, 15, 5};

    ASSERT_THAT(expected, ElementsAreArray(result));
    expected = {3, 3, 3, 3, 3};
    ASSERT_EQ(expected[init.rank()], part.getLocalSize());
}

TEST(CyclePartition, init_non_modulo){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    CyclePartition1D part(18, init);
    ASSERT_EQ(18, part.getGlobalSize());
    std::vector<unsigned int> result = part.getGlobalIndexRange();

    std::vector<unsigned int> expected(3);
    if (init.rank() == 0)
        expected = {0, 18, 5};
    else if (init.rank() == 1)
        expected = {1, 18, 5};
    else if (init.rank() == 2)
        expected = {2, 18, 5};
    else if (init.rank() == 3)
        expected = {3, 18, 5};
    else if (init.rank() == 4)
        expected = {4, 18, 5};

    ASSERT_THAT(expected, ElementsAreArray(result));
    expected = {4, 4, 4, 3, 3};
    ASSERT_EQ(expected[init.rank()], part.getLocalSize());
}

TEST(CyclePartition, local_to_global_index){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    CyclePartition1D part(18, init);
    std::vector<unsigned int> expected(part.getLocalSize());
    if (init.rank() == 0)
        expected = {0, 5, 10, 15};
    else if (init.rank() == 1)
        expected = {1, 6, 11, 16};
    else if (init.rank() == 2)
        expected = {2, 7, 12, 17};
    else if (init.rank() == 3)
        expected = {3, 8, 13};
    else if (init.rank() == 4)
        expected = {4, 9, 14};

    for (unsigned int local=0; local < part.getLocalSize(); local++)
        ASSERT_EQ(expected[local], part.localToGlobalIndex(local));
}

TEST(CyclePartition, global_to_local_index){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    CyclePartition1D part(18, init);
    std::vector<std::optional<unsigned int>> expected(part.getGlobalSize());
    if (init.rank() == 0)
        expected = {0, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 2, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 3, std::nullopt, std::nullopt };
    else if (init.rank() == 1)
        expected = {std::nullopt, 0, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 2, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 3, std::nullopt };
    else if (init.rank() == 2)
        expected = {std::nullopt, std::nullopt, 0, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 2, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 3 };
    else if (init.rank() == 3)
        expected = {std::nullopt, std::nullopt, std::nullopt, 0, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 2, std::nullopt, std::nullopt, std::nullopt, std::nullopt };
    else if (init.rank() == 4)
        expected = {std::nullopt, std::nullopt, std::nullopt, std::nullopt, 0, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 2, std::nullopt, std::nullopt, std::nullopt };

    for (unsigned int global=0; global < part.getGlobalSize(); global++)
        ASSERT_EQ(expected[global], part.globalToLocalIndex(global));

    // Example Usage for optional
    using ::testing::Optional;
    using ::testing::Eq;
    if (init.rank() == 1){
        ASSERT_THAT( part.globalToLocalIndex(6), Optional(1) );
        if (auto local = part.globalToLocalIndex(6)){
            ASSERT_EQ(*local, 1);
        }

        ASSERT_EQ( part.globalToLocalIndex(7), std::nullopt );
    }
}

TEST(CyclePartition, global_index_to_rank){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    CyclePartition1D part(18, init);
    std::vector<unsigned int> expected {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2};

    for (unsigned int global=0; global < part.getGlobalSize(); global++)
        ASSERT_EQ(expected[global], part.globalIndexToRank(global));
    
}

TEST(CyclePartition, partition_sizes){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    CyclePartition1D part(18, init);
    ASSERT_THAT(part.getPartitionSizes(), ElementsAreArray({4, 4, 4, 3, 3}));
}