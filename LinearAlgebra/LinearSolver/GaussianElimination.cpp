#include "GaussianElimination.h"

GaussianElimination::GaussianElimination(const LinearAlgebra &init){
    this->linalg = &init;
}

std::unique_ptr<AbstractVector> GaussianElimination::solve(const AbstractMatrix &A, const AbstractVector &b) const {
    std::unique_ptr<AbstractVector> result = A.mult(b);
    result->zeros();

    std::unique_ptr<AbstractMatrix> localA = A.clone();
    std::unique_ptr<AbstractVector> localb = b.clone();

    for (unsigned int col = 0; col < 1; col++){
        // header = Bcast( row corresponding to col)
        std::vector<double> header = localA->getRowValues(col);

        // for each local row (i)
        //      F = -1.0 * localA->getValue(i, col) / header[col]
        //      for each entry (j)
        //          localA->setLocalValue(i, j) = localA->getLocalValue(i,j) + F * header[j]

        unsigned int start_local = 0;


        std::vector<unsigned int> part = localA->getPartitionSize();
        for (unsigned int local_row = start_local; local_row < part[this->linalg->rank()]; local_row++){
            std::vector<double> localRow = localA->getLocalRowValues(local_row);
            double factor = -1.0 * localRow[col] / header[col];
            for (unsigned int entry = col; entry < localRow.size(); entry++){
                localRow[entry] += factor * header[entry];
            }
            localA->setLocalRowValues(local_row, localRow);
        }
    }

    localA->print();

    return result;
}