
#include "LinearAlgebra.h"

class Vector {
    public:
        Vector(unsigned int size, const LinearAlgebra& linalg);
        virtual ~Vector();

        void print() const;
        void add(const Vector &x);
        void setValues(const double &x);
        void scale(const double &x);

        unsigned int size() const {return this->global_size;}

        double getLocal(const unsigned int index) const;

        double length() const;
        
    private:
        unsigned int global_size, local_size;
        unsigned int global_starting_index;
        double* values;
        const LinearAlgebra* linalg;
};