#ifndef NUCLEUS_H_12LO
#define NUCLEUS_H_12LO

#include "gmock/gmock.h"
#include "mpi.h"

using ::testing::ElementsAreArray;

class Nucleus : public ::testing::EmptyTestEventListener {
    public:
        Nucleus();

        virtual void OnTestIterationStart(const ::testing::UnitTest& unit_test, int iteration) override;
        virtual void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration) override;
        virtual void OnEnvironmentsSetUpStart(const ::testing::UnitTest& unit_test) override;
        virtual void OnEnvironmentsTearDownStart(const ::testing::UnitTest& unit_test) override;
        virtual void OnTestCaseStart(const ::testing::TestCase& test_case) override;
        virtual void OnTestCaseEnd(const ::testing::TestCase& test_case) override;
        virtual void OnTestStart(const ::testing::TestInfo& info) override;
        virtual void OnTestEnd(const ::testing::TestInfo& info) override;
        virtual void OnTestPartResult(const ::testing::TestPartResult& result) override;

    private:
        std::string PrintFullTestCommentIfPresent(const ::testing::TestInfo& test_info);
        void OutputMessage(std::string message,
                           const unsigned int rank,
                           const unsigned int size,
                           const bool outputAllProcs = false) const;

        enum color {GREEN, RED, YELLOW};
        enum align {LEFT, RIGHT, CENTER, FULL};
        std::tuple<std::string, unsigned int, unsigned int> getMPIprefix(const color &c, const std::string &step, const align &alignment) const;
 
};

#endif // NUCLEUS_H_12LO
