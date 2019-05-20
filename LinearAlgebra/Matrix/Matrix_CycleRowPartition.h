#ifndef MATRIX_CYCLE_ROW_PARTITION_H_C98A
#define MATRIX_CYCLE_ROW_PARTITION_H_C98A

#include "AbstractMatrix.h"
#include "Vector_CyclePartition.h"

class Matrix_CycleRowPartition : public AbstractMatrix {

    public:
        Matrix_CycleRowPartition(unsigned int nRows,
            unsigned int nCols,
            const LinearAlgebra &linalg);
        Matrix_CycleRowPartition(const Matrix_CycleRowPartition &other);
        Matrix_CycleRowPartition& operator=(const Matrix_CycleRowPartition &other);
        virtual ~Matrix_CycleRowPartition();
        std::unique_ptr<AbstractMatrix> clone() const override;

        unsigned int nRows() const override { return this->nGlobalRows; }
        unsigned int nCols() const override { return this->nGlobalColumns; }
        std::vector<unsigned int> getPartitionSize() const override;
        IndexRange getGlobalRowIndexRange() const { return this->globalRowIndexRange; }
        std::vector<double> getRowValues(unsigned int row) const override;

        void setValues(const double &x) override;
        void setValues(const std::vector<std::vector<double>>& x) override;
        void setRowValues(const unsigned int row, const std::vector<double> &values) override;
        void zeros() override;

        unsigned int findRankWithIndex(const unsigned int index) const override;
        double getValue(const unsigned int row, const unsigned int col) const override;

        double frobeniusNorm() const override;

        std::unique_ptr<AbstractVector> mult(const AbstractVector &other) const override;

        void print(std::ostream &out = std::cout) const;

        void setLocalRowValues(const unsigned int row, const std::vector<double> &values) override;
        std::vector<double> getLocalRowValues(const unsigned int row) const override;

    private:
        unsigned int nGlobalRows, nGlobalColumns;
        unsigned int nLocalRows, nLocalColumns;
        const LinearAlgebra* linalg;
        IndexRange globalRowIndexRange;
        double** matrixStorage;
};

#endif // MATRIX_CYCLE_ROW_PARTITION_H_C98A
