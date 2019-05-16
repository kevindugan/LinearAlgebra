#ifndef NUCLEUS_H_12LO
#define NUCLEUS_H_12LO

#include "gmock/gmock.h"
#include "mpi.h"

using ::testing::ElementsAreArray;

class Nucleus : public ::testing::EmptyTestEventListener {
    public:
        Nucleus();

        virtual void OnTestStart(const ::testing::TestInfo& info) override;
        virtual void OnTestEnd(const ::testing::TestInfo& info) override;
        virtual void OnTestPartResult(const ::testing::TestPartResult& result) override;

    private:
        std::string PrintFullTestCommentIfPresent(const ::testing::TestInfo& test_info);
        void OutputMessage(std::string message, const unsigned int rank, const unsigned int size) const;
};

#endif // NUCLEUS_H_12LO
