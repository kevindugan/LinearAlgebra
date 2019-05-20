#ifndef VECTOR_CYCLE_PARTITION_H_C0A9
#define VECTOR_CYCLE_PARTITION_H_C0A9

#include "AbstractVector.h"

class Vector_CyclePartition : public AbstractVector {
    public:
        Vector_CyclePartition(unsigned int size, const LinearAlgebra& linalg);
        Vector_CyclePartition(const Vector_CyclePartition &other);
        Vector_CyclePartition& operator=(const Vector_CyclePartition &other);
        virtual ~Vector_CyclePartition();
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

        double getLocalValue(const unsigned int local_index) const override;
        IndexRange getGlobalIndexRange() const {return this->globalIndexRange;}
    
        void print(std::ostream &out = std::cout) const override;
        
    private:
        unsigned int local_size, global_size;
        double* values;
        const LinearAlgebra* linalg;
        IndexRange globalIndexRange;
};

#endif // VECTOR_CYCLE_PARTITION_H_C0A9
