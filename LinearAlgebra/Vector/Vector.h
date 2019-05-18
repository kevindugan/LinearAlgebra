#ifndef VECTOR_H_C99B
#define VECTOR_H_C99B
#include "LinearAlgebra.h"
#include <vector>
#include <utility>
#include "Partitioning.h"

class Vector {
    public:
        Vector(unsigned int size, const LinearAlgebra& linalg);
        Vector(const Vector &other);
        Vector& operator=(const Vector &other);
        virtual ~Vector();

        void print() const;
        void setValues(const double &x);
        void setValues(const std::vector<double> &x);
        void scale(const double &x);
        void zeros();

        Vector add(const double &scale, const Vector &other) const;
        Vector add(const Vector &other) const;
    
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

#endif // VECTOR_H_C99B
