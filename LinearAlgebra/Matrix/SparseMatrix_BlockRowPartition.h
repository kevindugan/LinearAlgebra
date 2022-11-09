#ifndef MATRIX_H_0C9S
#define MATRIX_H_0C9S

#include "AbstractMatrix.h"
#include "Vector_BlockPartition.h"

class SparseMatrix_BlockRowPartition : public AbstractMatrix {

    public:
        SparseMatrix_BlockRowPartition(unsigned int nRows,
               unsigned int nCols,
               const Parallel& linalg,
               const double tol=1.0e-16);
        SparseMatrix_BlockRowPartition(const SparseMatrix_BlockRowPartition &other);
        SparseMatrix_BlockRowPartition& operator=(const SparseMatrix_BlockRowPartition &other);
        virtual ~SparseMatrix_BlockRowPartition();
        std::unique_ptr<AbstractMatrix> clone() const override;

        unsigned int nRows() const override {return this->nGlobalRows;}
        unsigned int nCols() const override {return this->nGlobalColumns;}
        std::vector<unsigned int> getPartitionSize() const override;
        IndexRange getGlobalRowIndexRange() const override {return this->globalRowIndexRange;}
        std::vector<double> getRowValues(unsigned int row) const override;

        void setValues(const double& x) override;
        void setValues(const std::vector<std::vector<double>>& x) override;
        void setRowValues(const unsigned int row, const std::vector<double> &values) override;
        void zeros() override;

        unsigned int findRankWithIndex(const unsigned int index) const override;
        double getValue(const unsigned int row, const unsigned int col) const override;
        
        double frobeniusNorm() const override;
    
        std::unique_ptr<AbstractVector> mult(const AbstractVector &other) const override;

        void setLocalRowValues(const unsigned int row, const std::vector<double> &values) override;
        std::vector<double> getLocalRowValues(const unsigned int row) const override;

        void print(std::ostream &out = std::cout) const override;
        void printSparsity(std::ostream &out = std::cout) const;

    private:
        unsigned int nGlobalRows, nGlobalColumns;
        unsigned int nLocalRows, nLocalColumns;
        double zeroTol;
        const Parallel* linalg;
        IndexRange globalRowIndexRange;
        std::vector<double> values;
        std::vector<unsigned int> columnIndices, rowPtr;
};

#endif // MATRIX_H_0C9S
