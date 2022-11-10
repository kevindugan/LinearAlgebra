#include "Nucleus.h"
#include "BlockPartition1D.h"

TEST(BlockPartition, init_modulo){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    BlockPartition1D part(15, init);
    ASSERT_EQ(15, part.getGlobalSize());
    std::vector<unsigned int> result = part.getGlobalIndexRange();

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

    ASSERT_THAT(expected, ElementsAreArray(result));
    ASSERT_EQ(expected[1]-expected[0], part.getLocalSize());
}

TEST(BlockPartition, init_non_modulo){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    BlockPartition1D part(18, init);
    ASSERT_EQ(18, part.getGlobalSize());
    std::vector<unsigned int> result = part.getGlobalIndexRange();

    std::vector<unsigned int> expected(2);
    if (init.rank() == 0)
        expected = {0, 4};
    else if (init.rank() == 1)
        expected = {4, 8};
    else if (init.rank() == 2)
        expected = {8, 12};
    else if (init.rank() == 3)
        expected = {12, 15};
    else if (init.rank() == 4)
        expected = {15, 18};

    ASSERT_THAT(expected, ElementsAreArray(result));
    ASSERT_EQ(expected[1]-expected[0], part.getLocalSize());
}

TEST(BlockPartition, local_to_global_index){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    BlockPartition1D part(18, init);
    std::vector<unsigned int> expected(part.getLocalSize());
    if (init.rank() == 0)
        expected = {0, 1, 2, 3};
    else if (init.rank() == 1)
        expected = {4, 5, 6, 7};
    else if (init.rank() == 2)
        expected = {8, 9, 10, 11};
    else if (init.rank() == 3)
        expected = {12, 13, 14};
    else if (init.rank() == 4)
        expected = {15, 16, 17};

    for (unsigned int local=0; local < part.getLocalSize(); local++)
        ASSERT_EQ(expected[local], part.localToGlobalIndex(local));
}

TEST(BlockPartition, global_to_local_index){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    BlockPartition1D part(18, init);
    std::vector<std::optional<unsigned int>> expected(part.getGlobalSize());
    if (init.rank() == 0)
        expected = {0, 1, 2, 3, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt };
    else if (init.rank() == 1)
        expected = {std::nullopt, std::nullopt, std::nullopt, std::nullopt, 0, 1, 2, 3, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt };
    else if (init.rank() == 2)
        expected = {std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 0, 1, 2, 3, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt };
    else if (init.rank() == 3)
        expected = {std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 0, 1, 2, std::nullopt, std::nullopt, std::nullopt };
    else if (init.rank() == 4)
        expected = {std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 0, 1, 2 };

    for (unsigned int global=0; global < part.getGlobalSize(); global++)
        ASSERT_EQ(expected[global], part.globalToLocalIndex(global));

    // Example Usage for optional
    using ::testing::Optional;
    using ::testing::Eq;
    if (init.rank() == 1){
        ASSERT_THAT( part.globalToLocalIndex(5), Optional(1) );
        // Use this type of if block when referencing optional returns
        if (auto local = part.globalToLocalIndex(5)){
            ASSERT_EQ(*local, 1);
        }

        ASSERT_EQ( part.globalToLocalIndex(9), std::nullopt );
    }
}

TEST(BlockPartition, global_index_to_rank){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    BlockPartition1D part(18, init);
    std::vector<unsigned int> expected {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};

    for (unsigned int global=0; global < part.getGlobalSize(); global++)
        ASSERT_EQ(expected[global], part.globalIndexToRank(global));
    
}

TEST(BlockPartition, partition_sizes){
    Parallel init;
    ASSERT_EQ(init.size(), 5);

    BlockPartition1D part(18, init);
    ASSERT_THAT(part.getPartitionSizes(), ElementsAreArray({4, 4, 4, 3, 3}));
}