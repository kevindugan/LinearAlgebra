#ifndef ABSTRACT_MATRIX_H_QLWK
#define ABSTRACT_MATRIX_H_QLWK

#include "AbstractVector.h"

class AbstractMatrix {
    public:
        virtual ~AbstractMatrix() {}
        virtual std::unique_ptr<AbstractMatrix> clone() const = 0;

        virtual unsigned int nRows() const = 0;
        virtual unsigned int nCols() const = 0;
        virtual std::vector<unsigned int> getPartitionSize() const = 0;

        virtual void setValues(const double &x) = 0;
        virtual void setValues(const std::vector<std::vector<double>> &x) = 0;
        virtual void setRowValues(const unsigned int row, const std::vector<double> &values) = 0;
        virtual void zeros() = 0;

        virtual unsigned int findRankWithIndex(const unsigned int index) const = 0;
        virtual double getValue(const unsigned int row, const unsigned int col) const = 0;
        virtual std::vector<double> getRowValues(unsigned int row) const = 0;

        virtual double frobeniusNorm() const = 0;

        virtual std::unique_ptr<AbstractVector> mult(const AbstractVector &other) const = 0;

        virtual void print(std::ostream &out = std::cout) const {}

        virtual void setLocalRowValues(const unsigned int row, const std::vector<double> &values) = 0;
        virtual std::vector<double> getLocalRowValues(const unsigned int row) const = 0;
};

#endif // ABSTRACT_MATRIX_H_QLWK