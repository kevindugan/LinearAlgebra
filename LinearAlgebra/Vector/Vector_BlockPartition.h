#ifndef VECTOR_BLOCK_PARTITION_H_C99B
#define VECTOR_BLOCK_PARTITION_H_C99B

#include "AbstractVector.h"

class Vector_BlockPartition : public AbstractVector {
    public:
        Vector_BlockPartition(unsigned int size, const Parallel& linalg);
        Vector_BlockPartition(const Vector_BlockPartition &other);
        Vector_BlockPartition& operator=(const Vector_BlockPartition &other);
        virtual ~Vector_BlockPartition();
        std::unique_ptr<AbstractVector> clone() const;

        void setValues(const double &x) override;
        void setValues(const std::vector<double> &x) override;
        void scale(const double &x) override;
        void zeros() override;

        std::unique_ptr<AbstractVector> add(const AbstractVector &other) const override;
        std::unique_ptr<AbstractVector> add(const double &scale,
                                            const AbstractVector &other) const override;
    
        unsigned int size() const override { return this->global_size; }
        double l2norm() const override;

        double getValue(const unsigned int index) const override;
        void getGlobalValues(std::vector<double> &global_values) const override;

        std::vector<int> getPartitionSize() const override;
        unsigned int findRankWithIndex(const unsigned int index) const override;

        IndexRange getGlobalIndexRange() const {return this->globalIndexRange;}
        
        void print() const;

        double getLocalValue(const unsigned int local_index) const override;

    private:
        unsigned int local_size, global_size;
        double* values;
        const Parallel* linalg;
        IndexRange globalIndexRange;
};

#endif // VECTOR_BLOCK_PARTITION_H_C99B
