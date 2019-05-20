#ifndef ABSTRACT_VECTOR_H_C0Q9
#define ABSTRACT_VECTOR_H_C0Q9

#include "Partitioning.h"
#include "LinearAlgebra.h"
#include <vector>
#include <memory>

class AbstractVector {
    public:
        virtual ~AbstractVector() {}
        virtual std::unique_ptr<AbstractVector> clone() const = 0;

        virtual void setValues(const double &x) = 0;
        virtual void setValues(const std::vector<double> &x) = 0;
        virtual void scale(const double &x) = 0;
        virtual void zeros() = 0;

        virtual std::unique_ptr<AbstractVector> add(const double &scale, const AbstractVector &other) const = 0;
        virtual std::unique_ptr<AbstractVector> add(const AbstractVector &other) const = 0;

        virtual unsigned int size() const = 0;
        virtual double l2norm() const = 0;

        virtual double getValue(const unsigned int index) const = 0;
        virtual void getGlobalValues(std::vector<double> &global_values) const = 0;

        virtual std::vector<int> getPartitionSize() const = 0;
        virtual unsigned int findRankWithIndex(const unsigned int index) const = 0;

        virtual double getLocalValue(const unsigned int local_index) const = 0;
    
        virtual void print(std::ostream &out = std::cout) const {}
};

#endif // ABSTRACT_VECTOR_H_C0Q9
