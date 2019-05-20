#include "GaussianElimination.h"

GaussianElimination::GaussianElimination(){

}

std::unique_ptr<AbstractVector> GaussianElimination::solve(const AbstractMatrix &A, const AbstractVector &b) const {
    std::unique_ptr<AbstractVector> result = A.mult(b);
    result->zeros();

    std::unique_ptr<AbstractMatrix> localA = A.clone();
    std::unique_ptr<AbstractVector> localb = b.clone();

    return result;
}