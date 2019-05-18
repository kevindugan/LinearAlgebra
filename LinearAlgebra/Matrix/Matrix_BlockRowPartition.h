#ifndef MATRIX_H_0C9S
#define MATRIX_H_0C9S
#include "LinearAlgebra.h"
#include <vector>
#include "Partitioning.h"
#include "Vector_BlockPartition.h"

class Matrix_BlockRowPartition {

    public:
        Matrix_BlockRowPartition(unsigned int nRows,
               unsigned int nCols,
               const LinearAlgebra& linalg);
        virtual ~Matrix_BlockRowPartition();

        unsigned int nRows() const {return this->nGlobalRows;}
        unsigned int nCols() const {return this->nGlobalColumns;}
        std::vector<unsigned int> getPartitionSize() const;
        IndexRange getGlobalRowIndexRange() const {return this->globalRowIndexRange;}

        void setValues(const double& x);
        void setValues(const std::vector<std::vector<double>>& x);
        void zeros();
        double frobeniusNorm() const;
    
        std::unique_ptr<AbstractVector> mult(const AbstractVector &other) const;

    private:
        unsigned int nGlobalRows, nGlobalColumns;
        unsigned int nLocalRows, nLocalColumns;
        const LinearAlgebra* linalg;
        IndexRange globalRowIndexRange;
        double** matrixStorage;
};

#endif // MATRIX_H_0C9S
