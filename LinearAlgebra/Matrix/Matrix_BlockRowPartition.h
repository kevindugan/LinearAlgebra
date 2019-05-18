#ifndef MATRIX_H_0C9S
#define MATRIX_H_0C9S

#include "AbstractMatrix.h"
#include "Vector_BlockPartition.h"

class Matrix_BlockRowPartition : public AbstractMatrix {

    public:
        Matrix_BlockRowPartition(unsigned int nRows,
               unsigned int nCols,
               const LinearAlgebra& linalg);
        virtual ~Matrix_BlockRowPartition();

        unsigned int nRows() const override {return this->nGlobalRows;}
        unsigned int nCols() const override {return this->nGlobalColumns;}
        std::vector<unsigned int> getPartitionSize() const override;
        IndexRange getGlobalRowIndexRange() const {return this->globalRowIndexRange;}

        void setValues(const double& x) override;
        void setValues(const std::vector<std::vector<double>>& x) override;
        void zeros() override;
        
        double frobeniusNorm() const override;
    
        std::unique_ptr<AbstractVector> mult(const AbstractVector &other) const override;

    private:
        unsigned int nGlobalRows, nGlobalColumns;
        unsigned int nLocalRows, nLocalColumns;
        const LinearAlgebra* linalg;
        IndexRange globalRowIndexRange;
        double** matrixStorage;
};

#endif // MATRIX_H_0C9S
