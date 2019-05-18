#ifndef VECTOR_BLOCK_PARTITION_H_C99B
#define VECTOR_BLOCK_PARTITION_H_C99B
#include "LinearAlgebra.h"
#include <vector>
#include <utility>
#include "Partitioning.h"

class Vector_BlockPartition {
    public:
        Vector_BlockPartition(unsigned int size, const LinearAlgebra& linalg);
        Vector_BlockPartition(const Vector_BlockPartition &other);
        Vector_BlockPartition& operator=(const Vector_BlockPartition &other);
        virtual ~Vector_BlockPartition();

        void print() const;
        void setValues(const double &x);
        void setValues(const std::vector<double> &x);
        void scale(const double &x);
        void zeros();

        Vector_BlockPartition add(const double &scale, const Vector_BlockPartition &other) const;
        Vector_BlockPartition add(const Vector_BlockPartition &other) const;
    
        unsigned int size() const {return this->global_size;}

        unsigned int findRankWithIndex(const unsigned int index) const;
        double getValue(const unsigned int index) const;
        void getGlobalValues(std::vector<double> &global_values) const;

        double length() const;

        std::vector<int> getPartitionSize() const;

        IndexRange getGlobalIndexRange() const {return this->globalIndexRange;}
        
    private:
        unsigned int global_size, local_size;
        double* values;
        const LinearAlgebra* linalg;
        IndexRange globalIndexRange;
};

#endif // VECTOR_BLOCK_PARTITION_H_C99B
