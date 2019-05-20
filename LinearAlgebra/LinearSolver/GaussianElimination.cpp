#include "GaussianElimination.h"

GaussianElimination::GaussianElimination(){

}

std::unique_ptr<AbstractVector> GaussianElimination::solve(const AbstractMatrix &A, const AbstractVector &b) const {
    std::unique_ptr<AbstractVector> result = A.mult(b);
    result->zeros();

    std::unique_ptr<AbstractMatrix> localA = A.clone();
    std::unique_ptr<AbstractVector> localb = b.clone();

    for (unsigned int col = 0; col < localA->nCols(); col++){
        // header = Bcast( row corresponding to col)

        // for each local row (i)
        //      F = -1.0 * localA->getValue(i, col) / header[col]
        //      for each entry (j)
        //          localA->setLocalValue(i, j) = localA->getLocalValue(i,j) + F * header[j]
    }

    localA->print();

    return result;
}