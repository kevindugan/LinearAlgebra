#include "CyclePartition1D.h"
#include <math.h>

CyclePartition1D::CyclePartition1D(unsigned int size, const Parallel& comm)
{
    this->comm = &comm;
    this->global_size = size;
    Nucleus_ASSERT_GE(size, comm.size());

    // Calculate local size
    double ratio = float(this->global_size) / float(this->comm->size());
    unsigned int splitRank = this->global_size % this->comm->size();
    this->local_size = (this->comm->rank() < splitRank) ? ceil(ratio) : floor(ratio);

    // Calculate Global Index Range
    this->globalIndexRange.begin = this->comm->rank();
    this->globalIndexRange.end = this->global_size;
    this->globalIndexRange.skip = this->comm->size();
}

unsigned int CyclePartition1D::localToGlobalIndex(const unsigned int &local) const
{
    return local * this->globalIndexRange.skip + this->globalIndexRange.begin;
}

std::optional<unsigned int> CyclePartition1D::globalToLocalIndex(const unsigned int &global) const
{
    unsigned int local = (global - this->globalIndexRange.begin)/this->globalIndexRange.skip;
    return (this->comm->rank() == this->globalIndexToRank(global))
                ? std::optional<unsigned int>(local) : std::nullopt;
}

unsigned int CyclePartition1D::globalIndexToRank(const unsigned int &index) const
{
    Nucleus_ASSERT_LT(index, this->global_size)
    return index % this->comm->size();
}

std::vector<unsigned int> CyclePartition1D::getPartitionSizes() const
{
    std::vector<unsigned int> result(this->comm->size());
    MPI_Allgather(&this->local_size, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

std::vector<unsigned int> CyclePartition1D::getGlobalIndexRange() const
{
    return {this->globalIndexRange.begin, this->globalIndexRange.end, this->globalIndexRange.skip};
}