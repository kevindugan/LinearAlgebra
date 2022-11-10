#ifndef CYCLE_PARTITION_1D_H_5MQD
#define CYCLE_PARTITION_1D_H_5MQD

#include "Partition1D.h"
#include "Parallel.h"

class CyclePartition1D : public Partition1D {
    public:
        CyclePartition1D(unsigned int size, const Parallel& comm);

        unsigned int getGlobalSize() const { return this->global_size; }
        unsigned int getLocalSize() const { return this->local_size; }

        unsigned int localToGlobalIndex(const unsigned int &local) const;
        std::optional<unsigned int> globalToLocalIndex(const unsigned int &global) const;

        unsigned int globalIndexToRank(const unsigned int &index) const;
        std::vector<unsigned int> getPartitionSizes() const;
        std::vector<unsigned int> getGlobalIndexRange() const;

    private:
        unsigned int local_size, global_size;
        const Parallel* comm;
};

#endif // CYCLE_PARTITION_1D_H_5MQD