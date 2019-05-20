#ifndef LINEAR_SOLVER_H_09AS
#define LINEAR_SOLVER_H_09AS

#include "AbstractVector.h"
#include "AbstractMatrix.h"

class LinearSolver {

    public:
        virtual std::unique_ptr<AbstractVector> solve(const AbstractMatrix &A, const AbstractVector &b) const = 0;
};

#endif // LINEAR_SOLVER_H_09AS