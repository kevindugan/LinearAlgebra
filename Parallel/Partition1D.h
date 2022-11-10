#ifndef PARTITION_1D_H_309C
#define PARTITION_1D_H_309C

#include <vector>
#include <optional>

class Partition1D {
    public:
        virtual ~Partition1D() {}
        virtual unsigned int getGlobalSize() const = 0;
        virtual unsigned int getLocalSize() const = 0;

        virtual unsigned int localToGlobalIndex(const unsigned int &local) const = 0;
        virtual std::optional<unsigned int> globalToLocalIndex(const unsigned int &global) const = 0;

        virtual unsigned int globalIndexToRank(const unsigned int &index) const = 0;
        virtual std::vector<unsigned int> getPartitionSizes() const = 0;
        virtual std::vector<unsigned int> getGlobalIndexRange() const = 0;

    protected:
        // Index that is at the beginning of the range and and index one past the end
        struct IndexRange {
            IndexRange() : begin(0), end(1), skip(1) {}
            IndexRange(unsigned int a, unsigned int b) : begin(a), end(b), skip(1) {}
            IndexRange(unsigned int a, unsigned int b, unsigned int s) : begin(a), end(b), skip(s) {}
            unsigned int begin, end, skip;
        } globalIndexRange;
};

#endif // PARTITION_1D_H_309C