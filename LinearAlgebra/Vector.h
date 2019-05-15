
#include "LinearAlgebra.h"
#include <vector>
#include <utility>

// Index that is at the beginning of the range and and index one past the end
struct IndexRange {
    IndexRange() : begin(0), end(1) {}
    IndexRange(unsigned int a, unsigned int b) : begin(a), end(b) {}
    unsigned int begin, end;
};

class Vector {
    public:
        Vector(unsigned int size, const LinearAlgebra& linalg);
        virtual ~Vector();

        void print() const;
        void add(const Vector &x);
        void setValues(const double &x);
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