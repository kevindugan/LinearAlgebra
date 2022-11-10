#include "BlockPartition1D.h"
#include <math.h>

BlockPartition1D::BlockPartition1D(unsigned int size, const Parallel& comm)
{
    this->comm = &comm;
    this->global_size = size;
    Nucleus_ASSERT_GE(size, comm.size());

    // Calculate local size
    double ratio = float(this->global_size) / float(this->comm->size());
    unsigned int splitRank = this->global_size % this->comm->size();
    this->local_size = (this->comm->rank() < splitRank) ? ceil(ratio) : floor(ratio);

    // Calculate Global Index Range
    if (this->comm->rank() < splitRank){
        this->globalIndexRange.begin = this->comm->rank() * this->local_size;
        this->globalIndexRange.end = (this->comm->rank()+1) * this->local_size;
    } else {
        unsigned int offset = splitRank * (this->local_size + 1);
        this->globalIndexRange.begin = (this->comm->rank()-splitRank) * this->local_size + offset;
        this->globalIndexRange.end = (this->comm->rank()-splitRank +1) * this->local_size + offset;
    }
}

unsigned int BlockPartition1D::localToGlobalIndex(const unsigned int &local) const
{
    return local + this->globalIndexRange.begin;
}

std::optional<unsigned int> BlockPartition1D::globalToLocalIndex(const unsigned int &global) const
{
    unsigned int local = (global - this->globalIndexRange.begin)/this->globalIndexRange.skip;
    return (this->comm->rank() == this->globalIndexToRank(global))
                ? std::optional<unsigned int>(local) : std::nullopt;
}

unsigned int BlockPartition1D::globalIndexToRank(const unsigned int &index) const
{
    // Calculate index where change from s+1 to s sized partitions happens
    double ratio = float(this->global_size) / float(this->comm->size());
    unsigned int leftPerRank = ceil(ratio);
    unsigned int splitIndex = (this->global_size % this->comm->size()) * leftPerRank;
    // Calculate rank where index is in left split
    if (index < splitIndex)
        return floor( float(index)/float(leftPerRank) );
    // Or if index is in right split
    else {
        unsigned int rightPerRank = floor(ratio);
        return floor( float(index-splitIndex)/float(rightPerRank) ) + (this->global_size % this->comm->size());
    }
}

std::vector<unsigned int> BlockPartition1D::getPartitionSizes() const
{
    std::vector<unsigned int> result(this->comm->size());
    MPI_Allgather(&this->local_size, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

std::vector<unsigned int> BlockPartition1D::getGlobalIndexRange() const
{
    return {this->globalIndexRange.begin, this->globalIndexRange.end};
}