#ifndef VECTOR_H_C99B
#define VECTOR_H_C99B
#include "LinearAlgebra.h"
#include <vector>
#include <utility>
#include "Partitioning.h"

class Vector {
    public:
        Vector(unsigned int size, const LinearAlgebra& linalg);
        virtual ~Vector();

        void print() const;
        void add(const Vector &x);
        void setValues(const double &x);
        void setValues(const std::vector<double> &x);
        void scale(const double &x);
        void zeros();

        unsigned int size() const {return this->global_size;}

        double getLocal(const unsigned int index) const;

        double length() const;

        std::vector<unsigned int> getPartitionSize() const;

        IndexRange getGlobalIndexRange() const {return this->globalIndexRange;}
        
    private:
        unsigned int global_size, local_size;
        double* values;
        const LinearAlgebra* linalg;
        IndexRange globalIndexRange;
};

#endif // VECTOR_H_C99B