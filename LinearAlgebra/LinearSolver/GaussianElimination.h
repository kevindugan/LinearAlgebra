#ifndef GAUSSIAN_ELIMINATION_H_0C9A
#define GAUSSIAN_ELIMINATION_H_0C9A

#include "LinearSolver.h"

class GaussianElimination : public LinearSolver {

    public:
        GaussianElimination();

        std::unique_ptr<AbstractVector> solve(const AbstractMatrix &A, const AbstractVector &b) const override;
};

#endif // GAUSSIAN_ELIMINATION_H_0C9A