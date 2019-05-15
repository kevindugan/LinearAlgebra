
#include "LinearAlgebra.h"
#include <vector>

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
        
    private:
        unsigned int global_size, local_size;
        double* values;
        const LinearAlgebra* linalg;
};